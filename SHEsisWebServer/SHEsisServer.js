var sys=require('sys');
var express=require('express');
var app=express();
app.use(express.compress());
app.use(express.static('public'));

app.get('/',function(req,res){
console.log("request /");
res.redirect("SHEsis.html");
});

app.listen(8888);
console.log("SHEsis web server has started.");

