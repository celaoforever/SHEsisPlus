var shell=require('shelljs');
var sys=require('sys');
var express=require('express');
var app=express();
var bodyParser=require('body-parser');
var fs=require('fs');
var util=require('util');
var cp=require('child_process');
var path=require('path');
var Db=require('mongodb').Db;
var Server = require('mongodb').Server;
var kue = require('kue')
  , jobs = kue.createQueue()
  ;

app.use(express.compress());
app.use(express.static('public'));
app.use(bodyParser.json());
app.use(bodyParser.urlencoded());
bodyParser({strict:false});
app.use(bodyParser.json({limit: '50mb'}));
app.use(bodyParser.urlencoded({limit: '50mb', extended: true}));
app.get('/',function(req,res){
res.redirect("SHEsis.html");
var ip=req.connection.remoteAddress;
var time=getDateTimeFormated();
console.log(ip+" requested / at "+time);
insertToDB("get",{IP:ip,TIME:time});
});

app.post('/Results',function(req,res){
var ip=req.connection.remoteAddress;
var time=getDateTimeFormated();
console.log("post Results for "+ip+" at "+time);
Enqueue(req,res);
});

jobs.process('Run',12,function(job,done){
console.log("processing job",job.id);
var finalarg=job.data.arg.replace(/DuMmY1234_output/g,job.id);
RunSHEsis(finalarg,job.id,job.data.email,job.data.dataset);
done && done();
});

app.listen(5903);
console.log("SHEsis web server started at "+getDateTimeFormated());

function Enqueue(req,res){
var JSONARG={};
JSONARG.IP=req.connection.remoteAddress;
JSONARG.TIME=getDateTimeFormated();
var id=JSONARG.IP+"_"+getDateTime();
var args=setArgs(req,id,JSONARG);
insertToDB("post",JSONARG);
var job=jobs.create('Run',{arg:args,email:JSONARG.EMAIL,dataset:JSONARG.DATASET});
   job.on('complete',function(){
        var page="tmp/"+job.id+".html";
        setTimeout(function(){res.redirect(page);res.end();}, 1000);
        console.log('Job',job.id,'is done');
       
   })
      .on('failed',function(){
	console.log('Job',job.id,'has failed');
	})
   job.save();
}


function insertToDB(cl,content){
  var db=new Db("SHEsisLog",new Server('localhost',27017,{}),{w:-2,journal:false,fsync:false,safe: false});
db.open(function(err,db){
 var collection=db.collection(cl);
 collection.insert(content);
 db.close();
});
}

function setArgs(req,jobid,jsonarg){
var casedatafile="public/tmp/"+jobid+"_case.txt";
var ctrldatafile="public/tmp/"+jobid+"_ctrl.txt";
var qtldatafile="public/tmp/"+jobid+"_qtl.txt";
var output="public/tmp/DuMmY1234_output";
var args="";
var qtl=req.body.SelectPhenotype=="Case/Control"?false:true;
if(qtl){
	 fs.writeFile(qtldatafile,req.body.TextareaQTLdata,function(err){
             if(err){
                    console.log(err);
             };
        });
	args+=" --qtl --input "+qtldatafile;
	jsonarg.QTL=1;
}else{
	jsonarg.QTL=0;
	fs.writeFile(casedatafile,req.body.TextareaCasedata,function(err){
		if(err){
			console.log(err);
		};
	});

	fs.writeFile(ctrldatafile,req.body.TextareaControldata,function(err){
        	if(err){
                	console.log(err);
        	};
	});
	args+="--input-case "+casedatafile+" --input-ctrl "+ctrldatafile;
};
if(req.body.CheckBoxAnalysisTypeAssoc=="on"){
	args+=" --assoc";
	jsonarg.ASSOC=1;
}else
{
	jsonarg.ASSOC=0;
};
if(req.body.CheckBoxAnalysisTypeHWE=="on"){
	args+=" --hwe";
	jsonarg.HWE=1;
}else
{
	jsonarg.HWE=0;
};
if(req.body.CheckBoxAnalysisTypeHap=="on"){
	args+=" --haplo-EM";
	jsonarg.HAP=1;
}else
{
	jsonarg.HAP=0;
}
if(req.body.CheckBoxMultiCompPbased=="on"){
	args+=" --adjust";
	jsonarg.ADJUST=1;
}else
{
	jsonarg.ADJUST=0;
};


if(req.body.CheckBoxAnalysisTypeLD=="on"){
	if(req.body.SelectLDType=="Both case and control"){
		args+=" --ld";
		jsonarg.LD="LD in both";
	}else if(req.body.SelectLDType=="Just case"){
		args+=" --ld-in-case";
		jsonarg.LD="LD in case";
	}else if(req.body.SelectLDType=="Just control"){
		args+=" --ld-in-ctrl";
		jsonarg.LD="LD in ctrl";
	}else if(!req.body.SelectLDType){
		args+=" --ld";
		jsonarg.LD="LD in both";
	}else
	{
		jsonarg.LD=0;
	};
}
if(req.body.TextPermutation>0){
	args+=" --permutation "+req.body.TextPermutation;
};
jsonarg.PERMUTATION=req.body.TextPermutation;
args+=" --ploidy "+req.body.SelectPloidy;
args+=" --lft "+req.body.TextLFT;
args+=" --snpname-line \""+req.body.TextMarkername  + "\"";
args+=" --mask "+"\""+req.body.TextMask+"\"";
jsonarg.PLOIDY=req.body.SelectPloidy;
jsonarg.LFT=req.body.TextLFT;
jsonarg.SNP=req.body.TextMarkername;
jsonarg.MASK=req.body.TextMask;
args+=" --output "+output;
args+=" --webserver";
jsonarg.EMAIL=req.body.TextEmail;
jsonarg.DATASET=req.body.TextDatasetName;
return args;
};

app.get('/log', function(req, res) {
    console.log("request /log");
    var writeResp = function(err, result) {
        res.set('Content-Type', 'text/plain');
        if (err) {
            console.log("writeResp"+ err);
        } else {
            res.send(200, JSON.stringify(result));
            console.log("sent.");
        }
        res.end();
    };
    FindDB(req.param('table'), writeResp);
});

function FindDB(cl,callback){
  var db=new Db("SHEsisLog",new Server('localhost',27017,{}),{w:-2,journal:false,fsync:false,safe: false});
db.open(function(err,db){
if(err){console.log("error connecting to mongodb");}
console.log("connected to mongodb:"+cl);
 var collection=db.collection(cl);
        collection.find({}, {
            "_id": 0
        })
            .toArray(function(err, items) {
                callback(err, items);
                return;
            });
    }); 
}

function getDateTime() {

    var date = new Date();
    var hour = date.getHours();
    hour = (hour < 10 ? "0" : "") + hour;
    var min  = date.getMinutes();
    min = (min < 10 ? "0" : "") + min;
    var sec  = date.getSeconds();
    sec = (sec < 10 ? "0" : "") + sec;
    var year = date.getFullYear();
    var month = date.getMonth() + 1;
    month = (month < 10 ? "0" : "") + month;
    var day  = date.getDate();
    day = (day < 10 ? "0" : "") + day;
    return (year + month + day + hour +  min +  sec);

}


function getDateTimeFormated() {

    var date = new Date();
    var hour = date.getHours();
    hour = (hour < 10 ? "0" : "") + hour;
    var min  = date.getMinutes();
    min = (min < 10 ? "0" : "") + min;
    var sec  = date.getSeconds();
    sec = (sec < 10 ? "0" : "") + sec;
    var year = date.getFullYear();
    var month = date.getMonth() + 1;
    month = (month < 10 ? "0" : "") + month;
    var day  = date.getDate();
    day = (day < 10 ? "0" : "") + day;
    return (year+"-" + month +"-"+ day +" "+ hour +":"+  min +":"+  sec);

}
function SendEmail(addr,page,dataset){
    if(addr=="" || addr ==null)
	return;
console.log("sending email..");
var mail=
"Dear user, \n \
\n\
Your results are ready. You can access the results via "+page+" \n\
\n\
Please note that the results will be deleted after 5 days. If neccessary, please save the results before they become unavailable.\n\
\n\
Thank you for using SHEsis.\n\
\n\
Best Regards,\n\
SHEsis Team\n\
";
var cmd='echo "'+mail+'" '+'|mutt -s "SHEsis results';
if(dataset!=""){
	cmd+=" for dataset:"+dataset;
}
cmd+=' are ready" '+ addr +";";

var bin=cp.exec(cmd,function(err,stdout,stderr){
})
}

function SendInternalEmail(addr,sub,content){
if(addr=="" || addr == null)
	return;
console.log("Sending internal e-mail...");

var cmd= 'echo \''+content.replace(/\'/g,'')+'\' '+ '|mutt -s "'+sub+'" ' +addr+";";
//console.log(cmd)
var bin=cp.exec(cmd,function(err,stdout,stderr){
})
}


function RunSHEsis(args,jobid,email,dataset){
console.log('bin/SHEsis ',args);
var bin=cp.exec('bin/SHEsis '+args,{timeout:36000000},function(err,stdout,stderr){
var arrMatches=stdout.toString().match('ERROR.*');

if(err&&err.killed){
	    console.log("job",jobid,"timed out");
	    var msg="args: "+args;
	    SendInternalEmail("jiawei.shen@outlook.com","job "+jobid+ " timed out", msg);
}else if(arrMatches != null && arrMatches.length){
	    console.log("job",jobid,"failed");
	    var msg="args: "+args+"\n"+"msg: "+arrMatches[0];
	    SendInternalEmail("jiawei.shen@outlook.com","job "+jobid+" failed",msg);
}else if(err){	    
	    console.log("job",jobid,"error:",stderr);
	    var msg="args: "+args+"\n"+"SHEsis msg: "+stderr+"\nNode msg: "+err;
	    SendInternalEmail("jiawei.shen@outlook.com","job "+jobid+" error",msg); 
}else
{
     	SendEmail(email,"http://202.120.31.144:5903/tmp/"+jobid+".html",dataset);   
}
})
}
