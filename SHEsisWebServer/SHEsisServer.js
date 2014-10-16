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

app.get('/',function(req,res){
console.log("request /");
res.redirect("SHEsis.html");
var ip=req.connection.remoteAddress;
var time=getDateTimeFormated();
insertToDB("get",{IP:ip,TIME:time});
});

app.post('/Results',function(req,res){
console.log("post Results");
var ip=req.connection.remoteAddress;
var time=getDateTime();
var JSONARG={};
JSONARG.IP=ip;
JSONARG.TIME=getDateTimeFormated();
var args=setArgs(req,ip,time,JSONARG);
insertToDB("post",JSONARG);
console.log("cmd:bin/SHEsis "+args);
RunSHEsis(args,ip,time,res);
});

app.listen(5903);
console.log("SHEsis web server has started.");


function insertToDB(cl,content){
  var db=new Db("SHEsisLog",new Server('localhost',27017,{}),{w:-2,journal:false,fsync:false,safe: false});
db.open(function(err,db){
 var collection=db.collection(cl);
 collection.insert(content);
 db.close();
});
}

function setArgs(req,ip,time,jsonarg){
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

args+=" --ploidy "+req.body.SelectPloidy;
args+=" --lft "+req.body.TextLFT;
args+=" --snpname-line \""+req.body.TextMarkername  + "\"";
args+=" --mask "+"\""+req.body.TextMask+"\"";
jsonarg.PLOIDY=req.body.SelectPloidy;
jsonarg.LFT=req.body.TextLFT;
jsonarg.SNP=req.body.TextMarkername;
jsonarg.mask=req.body.TextMask;
args+=" --output "+output;
args+=" --webserver";
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


function RunSHEsis(args,ip,time,res){
var anyerr="";
 var bin=cp.exec('bin/SHEsis '+args,{timeout:3600000},function(err,stdout,stderr){
if(err){
if(err.killed){
 console.log("timed out");
 fs.readFile("public/Results.html",'utf8',function(err,data){
 if(err){
        return console.log(err);
 };
 results=data.replace(/SHOWRESULTSHERE/g,"<br><br><h2>Your request timed out. Please download the standalone version of SHEsis to run it on your local machine.</h2>");
 res.write(results);
 res.end();
})
}

else
 console.log("err");
return;
}
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
})
}


