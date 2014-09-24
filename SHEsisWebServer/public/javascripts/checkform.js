function validateAnalysisType(){
	var x=document.getElementsByName('assoc');
	if(x[0].checked){
	alert("Assoc checked");
	return false;
	};
};

function validateData(name){
	var x=document.getElementsByName(name);
	var a=x[0].value.split("\n");
	alert(a[0]);
	return false
}
function validateForm(){
validateAnalysisType();
validateData("casedata");
//validateData("controldata");
};
