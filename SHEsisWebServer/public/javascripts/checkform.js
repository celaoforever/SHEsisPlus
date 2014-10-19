function validateAnalysisType(){
	var count=0;
	if(document.getElementsByName('CheckBoxAnalysisTypeAssoc')[0].checked){
		count++;
	};
        if(document.getElementsByName('CheckBoxAnalysisTypeHWE')[0].checked){
                count++;
        };
        if(document.getElementsByName('CheckBoxAnalysisTypeHap')[0].checked){
                count++;
        };
        if(document.getElementsByName('CheckBoxAnalysisTypeLD')[0].checked){
                count++;
        };

	if(count == 0 ){
		setLegend("At least one analysis type should be checked.");
		return false;
	};
};

function setLegend(msg){
	var a=document.getElementsByName("legend");
	a[0].textContent= msg;
	a[0].style.color="red";
	a[0].style.fontWeight="bold";
	scroll(0,0);
}

function validateEmail(){
	var addr=document.getElementsByName("TextEmail")[0].value;
	if (addr == null || addr == "")
		return true;
	var atpos = addr.indexOf("@");
   	var dotpos = addr.lastIndexOf(".");
   	if (atpos< 1 || dotpos<atpos+2 || dotpos+2>=addr.length) {
      		setLegend("Invalid E-mail");  	
		return false;
	}
}

function validateData(name,exp){
        var ploidy=document.getElementsByName("SelectPloidy");
        var data=document.getElementsByName(name);
        var lines=data[0].value.split("\n");
        var expected=exp;
        var lineidx=0;
        for(i=0;i<lines.length;i++){

                var fileds=lines[i].trim().split(/[\s,]+/);
                if(fileds.length==1 && fileds[0]==""){
                        continue;
                }
                if(fileds.length==1){
                        setLegend("Error in "+(name=="TextareaCasedata"?"case data":"control data")+", line "+(i+1)+", Only sample ID found. No genotype data given.");
                        return false;
                };
                if(expected==0&&lineidx==0){
                        if((fileds.length-1)%ploidy[0].value != 0){
                                setLegend("Error in "+(name=="TextareaCasedata"?"case data":"control data")+", line "+ (i+1)+ ", Either the snp number or ploidy number is wrong.");
                                return false;
                        }else
                        {
                                expected=fileds.length;
                        };
                }else
                {
                        if(fileds.length!=expected)
                        {
                                setLegend("Error in "+(name=="TextareaCasedata"?"case data":"control data")+", line "+(i+1)+". Column number is not the same as that in the "+(exp==0?"first line.":"case data."));
                                return false;
                        };
                };
                lineidx++;
        };
        if(lineidx==0)
        {
                setLegend("Error in "+(name=="TextareaCasedata"?"case data":"control data")+". No data available.");
        return false;
        }
        return expected;

};

function validateQTLData(){
        var ploidy=document.getElementsByName("SelectPloidy");
        var data=document.getElementsByName("TextareaQTLdata");
        var lines=data[0].value.split("\n");
        var lineidx=0;
	var expected=0;
        for(i=0;i<lines.length;i++){
                var fileds=lines[i].trim().split(/[\s,]+/);
                if(fileds.length==1 && fileds[0]==""){
                        continue;
                }
                if(fileds.length<=2){
                        setLegend("Error in "+"line "+(i+1)+", no genotype data given.");
                        return false;
                };
                if(expected==0&&lineidx==0){
                        if((fileds.length-2)%ploidy[0].value != 0){
                                setLegend("Error in line "+ (i+1)+ ", Either the snp number or ploidy number is wrong.");
                                return false;
                        }else
                        {
                                expected=fileds.length;
                        };
                }else
                {
                        if(fileds.length!=expected)
                        {
                                setLegend("Error in line "+(i+1)+". Column number is not the same as that in the first line.");
                                return false;
                        };
                };
                lineidx++;
        };
        if(lineidx==0)
        {
                setLegend("Error: No data available.");
        return false;
        }
        return expected;
};
			 
function validateForm(){
//return ( validateAnalysisType()&&validateData("TextareaCasedata")&&validateData("TextareaControldata"));
//document.getElementsByName("SelectLDType")[0].disabled=false;

if(validateAnalysisType() == false){
return false;
};
if(validateEmail()==false){
return false;
}
if(document.getElementsByName("SelectPhenotype")[0].value=="Quantitative Trait"){
 if( validateQTLData() ==false){
	return false;
 }
}else{
var ret=validateData("TextareaCasedata",0);
if(ret==false){
return false;
};
if(validateData("TextareaControldata",ret)==false){
return false;
}
};
};
