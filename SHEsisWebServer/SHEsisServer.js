var sys=require('sys');
var express=require('express');
var app=express();
var bodyParser=require('body-parser');
var fs=require('fs');
var util=require('util');
var cp=require('child_process');
app.use(express.compress());
app.use(express.static('public'));
app.use(bodyParser.json());
app.use(bodyParser.urlencoded());



app.get('/',function(req,res){
console.log("request /");
res.redirect("SHEsis.html");
});

app.post('/Results',function(req,res){
console.log("post Results");
var ip=req.connection.remoteAddress;
var time=getDateTime();
var args=setArgs(req,ip,time);
//console.log(args);
var bin=cp.spawn('bin/SHEsis',[]);

var anyerr="";
bin.stdout.on("data",function(data){
	//var rePattern=new RegExp('/*ERROR.*$/');
	var arrMatches=data.toString().match('ERROR.*');
        if(arrMatches != null){
	if(arrMatches.length>0){
                anyerr=arrMatches[0];
		console.log(anyerr);
        };
	};
});
bin.on('exit',function(code){
    console.log('exited with code '+code);
});
fs.readFile("public/Results.html",'utf8',function(err,data){
 if(err){
	return console.log(err);
 }; 
var results="";
if(anyerr!=""){
 results=data.replace(/SHOWRESULTSHERE/g,anyerr);
}else
 {
  fs.readFile("tmp/"+ip+time+"output.html",'utf8',function(err,data){
 	if(err){
		return console.log(err);
	};
 results=data.replace(/SHOWRESULTSHERE/g,data); 
});
 }
 res.write(results);
 res.end();
});
//console.log(req.body);
});

function setArgs(req,ip,time){
var casedatafile="tmp/"+ip+time+"case.txt";
var ctrldatafile="tmp/"+ip+time+"ctrl.txt";
var output="tmp/"+ip+time+"output";
var args="";
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

if(req.body.CheckBoxAnalysisTypeAssoc=="on"){
	args+=" --assoc";
}
if(req.body.CheckBoxAnalysisTypeHWE=="on"){
	args+=" --hwe";
};
if(req.body.CheckBoxAnalysisTypeHap=="on"){
	args+=" --haplo";
};

if(req.body.CheckBoxAnalysisTypeLD=="on"){
	if(req.body.SelectLDType=="Both case and control"){
		args+=" --ld";
	}else if(req.body.SelectLDType=="Just case"){
		args+=" --ld-in-case";
	}else if(req.body.SelectLDType=="Just control"){
		args+=" --ld-in-ctrl";
	};
}

args+=" --ploidy "+req.body.SelectPloidy;
args+=" --lft "+req.body.TextLFT;
args+=" --snpname-line \""+req.body.TextMarkername  + "\"";
args+=" --mask "+"\""+req.body.TextMask+"\"";
args+=" --output "+output;
return args;
};

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
app.listen(8885);
console.log("SHEsis web server has started.");

