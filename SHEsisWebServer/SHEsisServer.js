var sys=require('sys');
var express=require('express');
var app=express();
var bodyParser=require('body-parser');
var fs=require('fs');
var util=require('util');
var cp=require('child_process');
var path=require('path');
app.use(express.compress());
app.use(express.static('public'));
app.use(bodyParser.json());
app.use(bodyParser.urlencoded());
bodyParser({strict:false});

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
var anyerr="";
console.log("cmd:bin/SHEsis "+args);
var bin=cp.exec('bin/SHEsis '+args,{timeout:6000000},function(err,stdout,stderr){
if(err){
 console.log(err);
}
//console.log(stdout);
var arrMatches=stdout.toString().match('ERROR.*');
if(arrMatches != null){
       if(arrMatches.length>0){
               anyerr=arrMatches[0];
               console.log(anyerr);
        };
};

var results="";
fs.readFile("public/Results.html",'utf8',function(err,data){
 if(err){
	return console.log(err);
 }; 
if(anyerr!=""){
 results=data.replace(/SHOWRESULTSHERE/g,anyerr);
 res.write(results);
 res.end();
}else
 {
 var filepath=ip+time+"output.html";
 
  fs.readFile(path.join("public/tmp/",filepath),'utf8',function(err,html){
 	if(err){
		return console.log(err);
	};
//})
 results=data.replace(/SHOWRESULTSHERE/g,html); 
  res.write(results);
 res.end();
});
 }
});
});
});

function setArgs(req,ip,time){
var casedatafile="public/tmp/"+ip+time+"case.txt";
var ctrldatafile="public/tmp/"+ip+time+"ctrl.txt";
var qtldatafile="public/tmp/"+ip+time+"qtl.txt";
var output="public/tmp/"+ip+time+"output";
var args="";
var qtl=req.body.SelectPhenotype=="Case/Control"?false:true;
if(qtl){
	 fs.writeFile(qtldatafile,req.body.TextareaQTLdata,function(err){
             if(err){
                    console.log(err);
             };
        });
	args+=" --qtl --input "+qtldatafile;
}else{

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
}
if(req.body.CheckBoxAnalysisTypeHWE=="on"){
	args+=" --hwe";
};
if(req.body.CheckBoxAnalysisTypeHap=="on"){
	args+=" --haplo-EM";
};

if(req.body.CheckBoxAnalysisTypeLD=="on"){
	if(req.body.SelectLDType=="Both case and control"){
		args+=" --ld";
	}else if(req.body.SelectLDType=="Just case"){
		args+=" --ld-in-case";
	}else if(req.body.SelectLDType=="Just control"){
		args+=" --ld-in-ctrl";
	}else if(!req.body.SelectLDType){
		args+=" --ld";
	};
}

args+=" --ploidy "+req.body.SelectPloidy;
args+=" --lft "+req.body.TextLFT;
args+=" --snpname-line \""+req.body.TextMarkername  + "\"";
args+=" --mask "+"\""+req.body.TextMask+"\"";
args+=" --output "+output;
args+=" --webserver";
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
app.listen(5903);
console.log("SHEsis web server has started.");

