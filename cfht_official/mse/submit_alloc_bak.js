/**
 * Attempt at using AJAX to "submit form and then plot SNR spectrum"
 * @author Nicolas FLAGEY, Billy Mahoney
 */

function isNumber(m)
{
    return typeof m == 'number' && !isNaN(m) && isFinite(m);
}


function updateQuerystring(submit_or_update)
{
    // a function to update the query-string when the user change a parameter
    
    error_qs = 0;
    
    /* Session ID */
    session_id  = document.getElementById("sessionID").value;
    querystring = "?sessionID="+session_id;
    
    /* Get spectro */
    mySpectro     = document.getElementById("spectro");
    mySpectroVal  = mySpectro.options[mySpectro.selectedIndex].value;
    /* Query string */
    querystring += "&spectro=" + mySpectroVal;

    /* Get targets parameters */
    myFile     = document.getElementById("upload");
    myFileVal  = myFile.files[0];
    myFovCtr     = document.getElementById("fov_ctr");
    myFovCtrVal  = myFovCtr.options[myFovCtr.selectedIndex].value;
    myFovCtrRA    = document.getElementById("fov_ctr_ra");
    myFovCtrRAVal = parseFloat(myFovCtrRA.value);
    myFovCtrDEC    = document.getElementById("fov_ctr_dec");
    myFovCtrDECVal = parseFloat(myFovCtrDEC.value);
    /* Query string */
    querystring += "&upload=" + myFileVal + "&fovctr=" + myFovCtrVal + "&fovctrra=" + myFovCtrRAVal + "&fovctrdec=" + myFovCtrDECVal;
    
    /* Get computation parameters */
    myMethod       = document.getElementById("meth");
    myMethodVal    = myMethod.options[myMethod.selectedIndex].value;
    myIterNum      = document.getElementById("fixiternum");
    myIterNumVal   = parseFloat(myIterNum.value);
    myAllocFrac    = document.getElementById("allocfrac");
    myAllocFracVal = parseFloat(myAllocFrac.value);
    /* Querystring */
    querystring += "&meth=" + myMethodVal + "&iternum=" + myIterNumVal + "&allocfrac=" + myAllocFracVal;
    
    
    /* Test values */
    if (!(isNumber(myIterNumVal)) || myIterNumVal<0 ) {
	alert("The number of iterations is invalid, please specify a positive number.");
	error_qs++;
    }
    if (!(isNumber(myAllocFracVal)) || myAllocFracVal<0  || myAllocFracVal>100 ) {
	alert("The allocation fraction is invalid, please specify a positive number smaller than 100.");
	error_qs++;
    }
    if (!(isNumber(myFovCtrRAVal)) || myFovCtrRAVal<0  || myFovCtrRAVal>360 ) {
	alert("The RA coordinate is invalid, please specify a positive number smaller than 360.");
	error_qs++;
    }
    if (!(isNumber(myFovCtrDECVal)) || myFovCtrDECVal<0  || myFovCtrDECVal>360 ) {
	alert("The DEC coordinate is invalid, please specify a positive number smaller than 360.");
	error_qs++;
    }
    if (!(isNumber(myAllocFracVal)) || myAllocFracVal<0  || myAllocFracVal>100 ) {
	alert("The allocation fraction is invalid, please specify a positive number smaller than 100.");
	error_qs++;
    }
    
}


function submit_alloc()
{
    
    /* Query all parameters (with 50ms delay for safety) */
    setTimeout(updateQuerystring(),50);

    /* If any error, stop here and disable Submit button, otherwise, go on */
    if (error_qs > 0) {
	document.getElementById("submit").disabled = true;
    }
    else {
	document.getElementById("submit").disabled = false;
	
	/* Change appearance of Submit button */
	document.getElementById("submit_color").bgColor = "#FF0000";
	document.getElementById("submit").value = 'Processing request ...';
	document.getElementById("submit").disabled = true;

	/* Prepare AJAX */
	var ajaxRequest;
	try {
	    ajaxRequest = new XMLHttpRequest();
	} catch(e) {
	    try {
		ajaxRequest = new ActiveXObject("Msxml2.XMLHTTP");
	    } catch (e) {
		try {
		    ajaxRequest = new ActiveXObject("Microsoft.XMLHTTP");
		} catch (e) {
		    alert ("Your browser cannot run this web application.");
		    return false;
		}
	    }
	}

	var form = document.getElementById('form_id');
	var formData = new FormData(form);
	/* ajaxRequest.open('POST', '../cgi-bin/mse/tmp/', true); */
	ajaxRequest.open('POST', "../cgi-bin/mse/mse_alloc_wrapper.py", true);
	ajaxRequest.send(formData);
	ajaxRequest.onreadystatechange = function() {
	    if (ajaxRequest.readyState == 4) {
		/* Printing results values */
		response = ajaxRequest.responseText;
		lines = response.split("\n");
		/* All text */
		alltext = "";
		for (i=0; i < lines.length; i++) {
		    alltext += lines[i] + "<br>"
		}
		cellParam = document.getElementById("results");
		cellParam.innerHTML = alltext;
	    } 
	}

	/* Call to mse_alloc_wrapper.py */
	ajaxRequest.open("GET","../cgi-bin/mse/mse_alloc_wrapper.py" + querystring, true);
	ajaxRequest.onreadystatechange = function() {
	    if (ajaxRequest.readyState == 4) {

		/* Printing results values */
		response = ajaxRequest.responseText;
		lines = response.split("\n");
		/* All text */
		alltext = "";
		for (i=1; i < lines.length; i++) {
		    alltext += lines[i] + "<br>"
		}
		
		cellParam = document.getElementById("results");
		cellParam.innerHTML = alltext;
		
		/* Restore appearance of Submit button */
		document.getElementById("submit_color").bgColor = "#FF0000";
		document.getElementById("submit").value = 'Submit';
		document.getElementById("submit").disabled = false;

	    } 
	}
	ajaxRequest.send(null);
    }
}



/* This is from https://www.html5rocks.com/en/tutorials/file/dndfiles/ */
function handleFileSelect(evt) {
    var files = evt.target.files; // FileList object
    var file = files[0]
    
    // files is a FileList of File objects. List some properties back into the HTML.
    var output = [];
    output.push('<li><strong>', escape(file.name), '</strong> (', file.type || 'n/a', ') - ',
                file.size, ' bytes, last modified: ',
                file.lastModifiedDate ? file.lastModifiedDate.toLocaleDateString() : 'n/a',
                '</li>');

    // Process only if ...
    if (!file.type.match('image.*')) {
	continue;
    }

    var reader = new FileReader();

    document.getElementById('tests').innerHTML = '<ul>' + output.join('') + '</ul>';
}

