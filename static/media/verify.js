function disable_upepsearch(thisform)
{
	thisform.min.disabled = true;	
	thisform.max.disabled = true;
	thisform.grace.disabled = true;
}

function enable_upepsearch(thisform)
{
	thisform.min.disabled = false;
	thisform.max.disabled = false;
	thisform.grace.disabled = false;
}

function toggle_upepsearch(thischeckbox){
if (thischeckbox.checked){
	enable_upepsearch(thischeckbox.form);
} else {
	disable_upepsearch(thischeckbox.form);
}
}

function add_option(thisform, thisselect, option, min, max)
{
if (!/^-?\d+$/.test(min.value)) {
	alert("The starting nucleotide of a marked region needs to be a valid positive integer.");
	min.focus();
	return false;
} 
if (!/^-?\d+$/.test(max.value)) {
	alert("The final nucleotide of a marked region needs to be a valid positive integer.");
	max.focus();
	return false;
}
var minnum = new Number(min.value);
var maxnum = new Number(max.value)
if (!(minnum < maxnum)){
	alert("The final nucleotide of a marked region needs to be after the starting nucleotide.");
	max.focus();
	return false;
}
if (!(minnum > 0)){
	alert("The starting nucleotide of a marked region needs to be a valid positive integer.");
	min.focus();
	return false;
}
if (/;/.test(option.value)){
	alert("Due to technical limitations, descriptions cannot contain \";\" characters.");
	option.focus();
	return false;
}  
var y=document.createElement('option');
y.text = "[" + min.value + ", " + max.value + "] " +option.value;
option.value = ""
min.value = ""
max.value = ""
try
  {
  thisselect.add(y,null); // standards compliant
  }
catch(ex)
  {
  thisselect.add(y); // IE only
  }
if (thisform.storage.value == ""){
	thisform.storage.value = y.text;
} else {
	thisform.storage.value += ("\n"+ y.text);
}
}

function remove_option(thisform, thisselect)
{
var temp=new Array();
temp = thisform.storage.value.split('\n');
while(thisselect.selectedIndex != -1){
	temp.splice(thisselect.selectedIndex, 1);
	thisselect.remove(thisselect.selectedIndex);
	}
thisform.storage.value="";

for (z in temp){
	thisform.storage.value += (temp[z] + "\n");
	}

if (!(thisform.storage.value == "")){
	thisform.storage.value = thisform.storage.value.slice(0,thisform.storage.value.length-2);
}
}

function matchsequence(sequence){
	for(i=0;i<sequence.length;i++){
		if (!(sequence.charAt(i)=='A'||sequence.charAt(i)=='T'||sequence.charAt(i)=='C'||sequence.charAt(i)=='G'||sequence.charAt(i)==' '||sequence.charAt(i)=='\n'||sequence.charAt(i)=='\r')){
		return false;
		}
	}
	return true;
}

function validate_conservedupepform(thisform)
{
thisform.seqquery.value=thisform.seqquery.value.toUpperCase();
if (thisform.seqquery.value==""){
	alert("Please enter a query sequence, or Genbank Accession number.");
	thisform.seqquery.focus();
	return false;
}
if (!/^-?\d+$/.test(thisform.match.value)) {
	animatedcollapse.show('uSettings');
	animatedcollapse.show('uAlignment');
	alert("'Alignment Options - Nucleotide Match' needs to be a valid integer.");
	thisform.match.focus();
	return false;
} 
if (!/^-?\d+$/.test(thisform.mismatch.value)) {
	animatedcollapse.show('uSettings');
	animatedcollapse.show('uAlignment');
	alert("'Alignment Options - Nucleotide Mismatch' needs to be a valid integer.");
	thisform.mismatch.focus();
	return false;
} 
if (!/^-?\d+$/.test(thisform.existence.value)) {
	animatedcollapse.show('uSettings');
	animatedcollapse.show('uAlignment');
	alert("'Alignment Options - Gap Existence Penalty' needs to be a valid integer.");
	thisform.existence.focus();
	return false;
} 
if (!/^-?\d+$/.test(thisform.extension.value)) {
	animatedcollapse.show('uSettings');
	animatedcollapse.show('uAlignment');
	alert("'Alignment Options - Gap Extension Penalty' needs to be a valid integer.");
	thisform.extension.focus();
	return false;
} 
if (!/^-?\d+$/.test(thisform.heatmapsize.value)) {
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hHeatmaps');
	alert("'Gradient Options - Heatmap Width' needs to be a valid positive integer.");
	thisform.heatmapsize.focus();
	return false;
} else {
	var varheatmapsize = new Number(thisform.heatmapsize.value);
	if (varheatmapsize>10000||varheatmapsize<400){
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Heatmap Width' needs to be between 400 and 10000");
		thisform.heatmapsize.focus();
		return false;
	}}
var count = 0;
var last = -1;
if (!thisform.gblack.value == "") {
	if (!/^-?\d+$/.test(thisform.gblack.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Black' needs to be a valid integer.");
		thisform.gblack.focus();
		return false;
	} else {
		var varblack = new Number(thisform.gblack.value);
		if (varblack>100||varblack<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Black' needs to be between 0 and 100");
			thisform.gblack.focus();
			return false;
		}
		last = varblack;
	}
	count++;
}
if (!thisform.ggreen.value == "") {
	if (!/^-?\d+$/.test(thisform.ggreen.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Green' needs to be a valid integer.");
		thisform.ggreen.focus();
		return false;
	} else {
		var vargreen = new Number(thisform.ggreen.value);
		if (vargreen>100||vargreen<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Green' needs to be between 0 and 100");
			thisform.ggreen.focus();
			return false;
		}
		if (!(vargreen > last)) {
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Green' must be greater than all preceding values");
			thisform.ggreen.focus();
			return false;
		}
		last = vargreen;
	}
	count++;
}
if (!thisform.gyellow.value == "") {
	if (!/^-?\d+$/.test(thisform.gyellow.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Yellow' needs to be a valid integer.");
		thisform.gyellow.focus();
		return false;
	} else {
		var varyellow = new Number(thisform.gyellow.value);
		if (varyellow>100||varyellow<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Yellow' needs to be between 0 and 100");
			thisform.gyellow.focus();
			return false;
		}
		if (!(varyellow > last)) {
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Yellow' must be greater than all preceding values");
			thisform.gyellow.focus();
			return false;
		}
		last = varyellow;	
	}
	count++;
}
if (!thisform.gorange.value == "") {
	if (!/^-?\d+$/.test(thisform.gorange.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Orange' needs to be a valid integer.");
		thisform.gorange.focus();
		return false;
	} else {
		var varorange = new Number(thisform.gorange.value);
		if (varorange>100||varorange<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Orange' needs to be between 0 and 100");
			thisform.gorange.focus();
			return false;
		}
		if (!(varorange > last)) {
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Orange' must be greater than all preceding values");
			thisform.gorange.focus();
			return false;
		}
		last = varorange;
	}
	count++;
}
if (!thisform.gred.value == "") {
	if (!/^-?\d+$/.test(thisform.gred.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Red' needs to be a valid integer.");
		thisform.gred.focus();
		return false;
	} else {
		var varred = new Number(thisform.gred.value);
		if (varred>100||varred<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Red' needs to be between 0 and 100");
			thisform.gred.focus();
			return false;
		}
		if (!(varred > last)) {
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Red' must be greater than all preceding values");
			thisform.gred.focus();
			return false;
		}
		last = varred;	
	}
	count++;
}
if (count < 2){
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hHeatmaps');
	alert("'Gradient Options' needs the values of at least two colours to create a gradient.");
	thisform.gblack.focus();
	return false;
}
if (!/^-?\d+$/.test(thisform.window.value)) {
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hHeatmaps');
	alert("'Window size' needs to be a valid positive integer.");
	thisform.window.focus();
	return false;
} else {
	var varwindow = new Number(thisform.window.value);
	if (!(varwindow>0)){
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Window size' needs to be a valid positive integer.");
		thisform.window.focus();
		return false;
	}
}
}


function validate_heatmapform(thisform)
{
thisform.alignquery.value=thisform.alignquery.value.toUpperCase();
thisform.alignref.value=thisform.alignref.value.toUpperCase();
if (thisform.alignquery.value==""){
	alert("Please enter a query sequence.");
	thisform.alignquery.focus();
	return false;
}
if (thisform.alignref.value==""){
	alert("Please enter a reference sequence.");
	thisform.alignref.focus();
	return false;
}
if (!(/^NM_\d+$/.test(thisform.alignquery.value)||/^XM_\d+$/.test(thisform.alignquery.value)||/^GI:\d+$/.test(thisform.alignquery.value))){
	if(!matchsequence(thisform.alignquery.value)){
	alert("The query sequence contains non-nucleotide characters. Please enter a valid query sequence.");
	thisform.alignquery.focus();
	return false;
	}
}
if (!(/^NM_\d+$/.test(thisform.alignref.value)||/^XM_\d+$/.test(thisform.alignref.value)||/^GI:\d+$/.test(thisform.alignref.value))){
	if(!matchsequence(thisform.alignref.value)){
	alert("The reference sequence contains non-nucleotide characters. Please enter a valid reference sequence.");
	thisform.alignref.focus();
	return false;
	}
}
if (!/^-?\d+$/.test(thisform.match.value)) {
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hAlignment');
	alert("'Alignment Options - Nucleotide Match' needs to be a valid integer.");
	thisform.match.focus();
	return false;
} 
if (!/^-?\d+$/.test(thisform.mismatch.value)) {
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hAlignment');
	alert("'Alignment Options - Nucleotide Mismatch' needs to be a valid integer.");
	thisform.mismatch.focus();
	return false;
} 
if (!/^-?\d+$/.test(thisform.existence.value)) {
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hAlignment');
	alert("'Alignment Options - Gap Existence Penalty' needs to be a valid integer.");
	thisform.existence.focus();
	return false;
} 
if (!/^-?\d+$/.test(thisform.extension.value)) {
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hAlignment');
	alert("'Alignment Options - Gap Extension Penalty' needs to be a valid integer.");
	thisform.extension.focus();
	return false;
}
if (!/^-?\d+$/.test(thisform.heatmapsize.value)) {
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hHeatmaps');
	alert("'Gradient Options - Heatmap Width' needs to be a valid positive integer.");
	thisform.heatmapsize.focus();
	return false;
} else {
	var varheatmapsize = new Number(thisform.heatmapsize.value);
	if (varheatmapsize>10000||varheatmapsize<400){
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Heatmap Width' needs to be between 400 and 10000");
		thisform.heatmapsize.focus();
		return false;
	}}
var count = 0;
var last = -1;
if (!thisform.gblack.value == "") {
	if (!/^-?\d+$/.test(thisform.gblack.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Black' needs to be a valid integer.");
		thisform.gblack.focus();
		return false;
	} else {
		var varblack = new Number(thisform.gblack.value);
		if (varblack>100||varblack<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Black' needs to be between 0 and 100");
			thisform.gblack.focus();
			return false;
		}
		last = varblack;
	}
	count++;
}
if (!thisform.ggreen.value == "") {
	if (!/^-?\d+$/.test(thisform.ggreen.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Green' needs to be a valid integer.");
		thisform.ggreen.focus();
		return false;
	} else {
		var vargreen = new Number(thisform.ggreen.value);
		if (vargreen>100||vargreen<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Green' needs to be between 0 and 100");
			thisform.ggreen.focus();
			return false;
		}
		if (!(vargreen > last)) {
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Green' must be greater than all preceding values");
			thisform.ggreen.focus();
			return false;
		}
		last = vargreen;
	}
	count++;
}
if (!thisform.gyellow.value == "") {
	if (!/^-?\d+$/.test(thisform.gyellow.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Yellow' needs to be a valid integer.");
		thisform.gyellow.focus();
		return false;
	} else {
		var varyellow = new Number(thisform.gyellow.value);
		if (varyellow>100||varyellow<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Yellow' needs to be between 0 and 100");
			thisform.gyellow.focus();
			return false;
		}
		if (!(varyellow > last)) {
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Yellow' must be greater than all preceding values");
			thisform.gyellow.focus();
			return false;
		}
		last = varyellow;	
	}
	count++;
}
if (!thisform.gorange.value == "") {
	if (!/^-?\d+$/.test(thisform.gorange.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Orange' needs to be a valid integer.");
		thisform.gorange.focus();
		return false;
	} else {
		var varorange = new Number(thisform.gorange.value);
		if (varorange>100||varorange<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Orange' needs to be between 0 and 100");
			thisform.gorange.focus();
			return false;
		}
		if (!(varorange > last)) {
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Orange' must be greater than all preceding values");
			thisform.gorange.focus();
			return false;
		}
		last = varorange;
	}
	count++;
}
if (!thisform.gred.value == "") {
	if (!/^-?\d+$/.test(thisform.gred.value)) {
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Gradient Options - Red' needs to be a valid integer.");
		thisform.gred.focus();
		return false;
	} else {
		var varred = new Number(thisform.gred.value);
		if (varred>100||varred<0){
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Red' needs to be between 0 and 100");
			thisform.gred.focus();
			return false;
		}
		if (!(varred > last)) {
			animatedcollapse.show('hSettings');
			animatedcollapse.show('hHeatmaps');
			alert("'Gradient Options - Red' must be greater than all preceding values");
			thisform.gred.focus();
			return false;
		}
		last = varred;	
	}
	count++;
}
if (count < 2){
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hHeatmaps');
	alert("'Gradient Options' needs the values of at least two colours to create a gradient.");
	thisform.gblack.focus();
	return false;
}
if (!/^-?\d+$/.test(thisform.window.value)) {
	animatedcollapse.show('hSettings');
	animatedcollapse.show('hHeatmaps');
	alert("'Window size' needs to be a valid positive integer.");
	thisform.window.focus();
	return false;
} else {
	var varwindow = new Number(thisform.window.value);
	if (!(varwindow>0)){
		animatedcollapse.show('hSettings');
		animatedcollapse.show('hHeatmaps');
		alert("'Window size' needs to be a valid positive integer.");
		thisform.window.focus();
		return false;
	}
}
for (i=0;i<thisform.markers.length;i++) {
    	thisform.markers.options[i].selected = true;
}
}

function reset_heatmapform(thisform){
	var doit = confirm('This will reset this page. Continue?');
	if (doit){
		for (i=0;i<thisform.markers.length;i++) {
    			thisform.markers.options[i].selected = true;
		}
		remove_option(thisform, thisform.markers);		
		return true;
	} else {return false;}
}

function quick_add(thisselect, value){
	y=document.createElement('option');
	y.text = value;
	try {
 	    thisselect.add(y,null); // standards compliant
	}
	catch(ex) {
            thisselect.add(y); // IE only
	}
}

function init(){
thisform = document.getElementById('heatmapform')
if (!(thisform.storage.value == "")){
	temp = thisform.storage.value.split('\n');
	for (x in temp){
		quick_add(thisform.markers, temp[x]);
	}
}
}


window.onload=init;
