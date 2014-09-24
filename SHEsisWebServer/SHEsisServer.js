var sys=require('sys');
var express=require('express');
var app=express();
var bodyParser=require('body-parser');
var fs=require('fs');
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
fs.readFile("public/Results.html",'utf8',function(err,data){
 if(err){
	return console.log(err);
 };
 var results=data.replace(/SHOWRESULTSHERE/g,req.body.controldata);
 res.write(results);
 res.end();
});
console.log(req.body.controldata);
});


app.listen(8888);
console.log("SHEsis web server has started.");

