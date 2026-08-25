import re
import os

# =======================================================
# 1. Patch HaplotypeLD.cpp (Fix mask indexing logic)
# =======================================================
with open("HaplotypeLD.cpp", "r") as f:
    c = f.read()

new_phaseAll = """void HaplotypeLD::phaseAll(){
	int numSnp = this->SnpIdx.size();
	int startsnp, endsnp;
	int count = 0;
	int curSnpNum = numSnp;
	while(true)
	{
		int numOfBlock = (curSnpNum - 1) / this->numSnpPerBlock + 1;
		this->PhasedData.reset(new SHEsisData(this->data->getSampleNum(), numOfBlock, this->data->getNumOfChrSet()));
		this->PhasedData->vLabel = this->data->vLabel;
		this->PhasedData->vQuantitativeTrait = this->data->vQuantitativeTrait;
		for(int i = 0; i < numOfBlock; i++){
			startsnp = i * this->numSnpPerBlock;
			endsnp = MIN((1 + i) * this->numSnpPerBlock - 1, curSnpNum - 1);
			if (count == 0) {
				this->getHaplotypeSub(startsnp, endsnp, this->data);
			} else {
				this->getHaplotypeSub(startsnp, endsnp, this->DataToPhase);
			}
			this->reducePhasedHap(count % 2 ? this->reduced : this->hapcode, this->PhasedData, this->hap->PhasedData, i);
		}
		this->curCode = 1;
		if (!this->hapcode.empty() && !this->reduced.empty()){
			this->updateHapCode(count % 2 ? this->hapcode : this->reduced, count % 2 ? this->reduced : this->hapcode);
		}
		if (1 == numOfBlock){
			this->hapcode = count % 2 ? this->reduced : this->hapcode;
			break;
		}
		this->DataToPhase.reset(new SHEsisData(this->data->getSampleNum(), numOfBlock, this->data->getNumOfChrSet()));
		this->DataToPhase->vLabel = this->data->vLabel;
		this->DataToPhase->vQuantitativeTrait = this->data->vQuantitativeTrait;
		this->DataToPhase->mGenotype.resize(extension(this->PhasedData->mGenotype));
		this->DataToPhase->mGenotype = this->PhasedData->mGenotype;
		curSnpNum = numOfBlock;
		count++;
	};
	this->hap.reset(new HaplotypeEM(this->data));
	this->hap->PhasedData.resize(extension(this->PhasedData->mGenotype));
	this->hap->PhasedData=this->PhasedData->mGenotype;
	this->hap->SnpIdx.clear();
	this->hap->setFreqThreshold(this->lft);
	this->hap->SnpIdx.push_back(0);
	this->hap->getResults();
	for(int i=0;i<this->hap->Results.haplotypes.size();i++){
		std::string ss;
		boost::shared_ptr<short[]> cur=this->hap->Results.haplotypes[i];
		boost::shared_ptr<short[]> res(new short[this->SnpIdx.size()]);
		int idx=0;
		for(int j=0;j<1;j++){
			hapMap::map_by<code>::const_iterator iter=this->hapcode.by<code>().find(cur[j]);
			BOOST_ASSERT(iter!=this->hapcode.by<code>().end());
			ss=iter->get<haplo>();
			std::vector<std::string> vec;
			boost::algorithm::split(vec,ss,boost::algorithm::is_any_of(","));
			for(int k=0;k<vec.size()-1;k++){
				res[idx++]=(short)atoi(vec[k].data());
			}
		}
		this->haplotypes.push_back(res);
	}
}"""

new_getSub = """void HaplotypeLD::getHaplotypeSub(int start, int end, boost::shared_ptr<SHEsisData> d){
	std::vector<short> m(d->getSnpNum(), 0);
	if (d == this->data) {
		for(int i = start; i <= end; i++){
			m[this->SnpIdx[i]] = 1;
		};
	} else {
		for(int i = start; i <= end; i++){
			m[i] = 1;
		};
	}

	int maxIteration = 50;
	int minIteraction = 15;
	int iteration = 0;
	std::vector<boost::shared_ptr<HaplotypeEM> > locals;
	std::vector<long> checksums;
	std::vector<int> numOfHaps;
	int finalchoose = -1;
	while(iteration < maxIteration){
		boost::shared_ptr<HaplotypeEM> local(new HaplotypeEM(d, end - start + 1, m));
		local->setSeed(iteration);
		local->setFreqThreshold(this->lft);
		local->startHaplotypeAnalysis();
		locals.push_back(local);
		checksums.push_back(getCheckSum(local->PhasedData));
		numOfHaps.push_back(local->Results.haplotypes.size());
		iteration++;
		if(iteration > minIteraction){
			finalchoose = getSelection(checksums, numOfHaps);
			if(-1 != finalchoose)
				break;
		}
	}
	if(-1 == finalchoose)
		finalchoose = getSelection(checksums, numOfHaps, false);
	this->hap = locals[finalchoose];
}"""

c = re.sub(r'void\s+(?:SHEsis::)?HaplotypeLD::phaseAll\s*\(\)\s*\{[\s\S]*?\n\}', new_phaseAll, c)
c = re.sub(r'void\s+(?:SHEsis::)?HaplotypeLD::getHaplotypeSub\s*\(int\s+start,\s*int\s+end,\s*boost::shared_ptr<SHEsisData>\s+d\)\s*\{[\s\S]*?\n\}', new_getSub, c)

with open("HaplotypeLD.cpp", "w") as f:
    f.write(c)

# =======================================================
# 2. Patch logistic.cpp (mlpack 2.x API & p-value call)
# =======================================================
with open("logistic.cpp", "r") as f:
    c = f.read()

new_log = """void logistic::regress()
{
    BOOST_ASSERT(this->regressors.n_cols > 0 && this->regressors.n_rows > 0 && this->responses.size() > 0);
    BOOST_ASSERT(this->regressors.n_cols == this->responses.size());
    arma::Row<size_t> uresp = arma::conv_to<arma::Row<size_t>>::from(this->responses.t());
    mlpack::regression::LogisticRegressionFunction<arma::mat> lrf(this->regressors, uresp, this->lambda);
    if (this->optimizerType == "sgd") {
        mlpack::optimization::SGD<mlpack::regression::LogisticRegressionFunction<arma::mat>> sgdOpt(lrf);
        sgdOpt.MaxIterations() = this->maxIterations;
        sgdOpt.Tolerance() = this->tolerance;
        sgdOpt.StepSize() = 0.01;
        mlpack::regression::LogisticRegression<arma::mat> lr(sgdOpt);
        this->coef = lr.Parameters();
    } else {
        mlpack::optimization::L_BFGS<mlpack::regression::LogisticRegressionFunction<arma::mat>> lbfgsOpt(lrf);
        lbfgsOpt.MaxIterations() = this->maxIterations;
        lbfgsOpt.MinGradientNorm() = this->tolerance;
        mlpack::regression::LogisticRegression<arma::mat> lr(lbfgsOpt);
        this->coef = lr.Parameters();
    }
    this->getPvalue();
}"""

c = re.sub(r'void\s+(?:SHEsis::)?logistic::regress\s*\(\)\s*\{[\s\S]*?\n\}', new_log, c)
with open("logistic.cpp", "w") as f:
    f.write(c)

# =======================================================
# 3. Patch linear.cpp (mlpack 2.x API & p-value call)
# =======================================================
with open("linear.cpp", "r") as f:
    c = f.read()

new_lin = """void linear::regress()
{
    BOOST_ASSERT(this->regressors.n_cols > 0 && this->regressors.n_rows > 0 && this->responses.size() > 0);
    BOOST_ASSERT(this->regressors.n_cols == this->responses.size());
    mlpack::regression::LinearRegression lr(this->regressors, this->responses, this->lambda);
    this->coef = lr.Parameters();
    lr.Predict(this->regressors, this->predictions);
    this->getPvalue();
}"""

c = re.sub(r'void\s+(?:SHEsis::)?linear::regress\s*\(\)\s*\{[\s\S]*?\n\}', new_lin, c)
with open("linear.cpp", "w") as f:
    f.write(c)

# =======================================================
# 4. Patch SHEsisWebServer/SHEsisServer.js (Docker network)
# =======================================================
server_path = "SHEsisWebServer/SHEsisServer.js"
if os.path.exists(server_path):
    with open(server_path, "r") as f:
        code = f.read()

    kue_patch = """
var redisHost = process.env.REDIS_HOST || "127.0.0.1";
var redisPort = parseInt(process.env.REDIS_PORT || 6379, 10);
var kue = require("kue");
kue.redis.createClient = function() {
  return require("redis").createClient(redisPort, redisHost);
};
var jobs = kue.createQueue({
  redis: { host: redisHost, port: redisPort }
});
"""
    code = re.sub(r'var\s+kue\s*=\s*require\([^\)]+\)[\s\S]*?jobs\s*=\s*kue\.createQueue\(\)[\s\S]*?;', kue_patch, code)
    code = code.replace("'localhost',27017", "process.env.MONGO_HOST || 'localhost', parseInt(process.env.MONGO_PORT || 27017, 10)")
    code = code.replace('"localhost",27017', "process.env.MONGO_HOST || 'localhost', parseInt(process.env.MONGO_PORT || 27017, 10)")
    code = code.replace("w:-2", "w:1")
    code = re.sub(r'fs\.writeFile\((datafile|covarfile),\s*([^\s,]+),\s*function\([^)]*\)\s*\{[\s\S]*?\}\);', r'fs.writeFileSync(\1, \2 || "");', code)
    code = code.replace('var id=JSONARG.IP+"_"+getDateTime();', 'var id=JSONARG.IP.replace(/[^0-9a-zA-Z]/g, "_")+"_"+getDateTime();')

    with open(server_path, "w") as f:
        f.write(code)

print("SHEsisPlus patch applied successfully.")
