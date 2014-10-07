function ResetForm(){
document.getElementsByName("TableInputData")[0].style.display="";
document.getElementsByName("TableInputDataQTL")[0].style.display="none";
}


function selectQTL(){
 var elementPheno = document.getElementsByName("SelectPhenotype");
 var elementLD    = document.getElementsByName("SelectLDType");
 var elementInputData=document.getElementsByName("TableInputData");
 var elementInputDataQTL=document.getElementsByName("TableInputDataQTL");
 if(elementPheno[0].value == "Quantitative Trait"){
 	elementLD[0].value="Both case and control";
 	elementLD[0].disabled=true;
	elementInputData[0].style.display="none";
	elementInputDataQTL[0].style.display="";
 }else{
	elementLD[0].disabled=false;
	elementInputData[0].style.display="";
	elementInputDataQTL[0].style.display="none";
}
}

function loadSampleBinary(){
var casedata="id1 1 2 A G G C C T\nid2 1 1 A A 0 0 C C\nid3 1 1 A A G G C T\nid4 2 1 G G G G T C\nid5 2 2 A G G C T T\nid6 3 3 A G G C C C\nid7 3 2 A A G G C C\nid8 0 0 G G C G C C\nid9 2 2 A A C G T T";
var ctrldata="id11 2 2 A A G C C T\nid12 1 1 A A C C C C\nid13 1 1 A A G G C T\nid14 1 1 G G G G T T\nid15 1 2 A A G C T T \nid16 1 3 A A G C T C \nid17 3 2 A A G G C C\nid18 2 1 G G C G C T\nid19 2 2 A A C G 0 0";
document.getElementsByName("TextareaCasedata")[0].value=casedata;
document.getElementsByName("TextareaControldata")[0].value=ctrldata;
document.getElementsByName('CheckBoxAnalysisTypeAssoc')[0].checked=true;
document.getElementsByName('CheckBoxAnalysisTypeHWE')[0].checked=true;
document.getElementsByName('CheckBoxAnalysisTypeHap')[0].checked=true;
document.getElementsByName('CheckBoxAnalysisTypeLD')[0].checked=true;
document.getElementsByName("SelectPhenotype")[0].value="Case/Control";
document.getElementsByName("SelectPloidy")[0].value="2";
document.getElementsByName("SelectLDType")[0].value="Both case and control";
document.getElementsByName("TextMarkername")[0].value="rs11,rs22,rs33,rs44";
document.getElementsByName("TextLFT")[0].value="0.03";
document.getElementsByName("TextMask")[0].value="1,1,0,1";
};

function loadSampleQTL(){
var qtldata="id1 2.3 1 2 A G G C C T\nid2 1.3 1 1 A A 0 0 C C\nid3 2.1 1 1 A A G G C T\nid4 3.2 2 1 G G G G T C\nid5 3.45 2 2 A G G C T T\nid6 3.33 3 3 A G G C C C\nid7 5.43 3 2 A A G G C C\nid8 3.21 0 0 G G C G C C\nid9 2.22 2 2 A A C G T T\nid10 3.2 1 1 A A C C T T\nid11 3.22 2 2 A A G C C T\nid12 4.56 1 1 A A C C C C\nid13 5.4 1 1 A A G G C T\nid14 2.24 1 1 G G G G T T\nid15 4.55 1 2 A A G C T T \nid16 6.5 1 3 A A G C T C \nid17 7.6 3 2 A A G G C C\nid18 4.6 2 1 G G C G C T\nid19 5.5 2 2 A A C G 0 0";
document.getElementsByName("TextareaQTLdata")[0].value=qtldata;
document.getElementsByName('CheckBoxAnalysisTypeAssoc')[0].checked=true;
document.getElementsByName('CheckBoxAnalysisTypeHWE')[0].checked=true;
document.getElementsByName('CheckBoxAnalysisTypeHap')[0].checked=true;
document.getElementsByName('CheckBoxAnalysisTypeLD')[0].checked=true;
document.getElementsByName("SelectPhenotype")[0].value="Quantitative Trait";
document.getElementsByName("SelectPloidy")[0].value="2";
document.getElementsByName("SelectLDType")[0].value="Both case and control";
document.getElementsByName("TextMarkername")[0].value="rs11,rs22,rs33,rs44";
document.getElementsByName("TextLFT")[0].value="0.03";
document.getElementsByName("TextMask")[0].value="1,1,0,1";

}
