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
    
    /* Get instrument parameters */
    /*   telescope */
    /* myCoating     = document.getElementById("coating");
    myCoatingVal  = myCoating.options[myCoating.selectedIndex].value; */
    /* Query string */
    /* querystring += "&coating=" + myCoatingVal; */
    querystring += "&coating=ZeCoat";
    
    
    /* Get observing conditions */
    mySeeing     = document.getElementById("seeing");
    mySeeingVal  = parseFloat(mySeeing.value);
    myAirmass    = document.getElementById("airmass");
    myAirmassVal = parseFloat(myAirmass.value);
    mySkyMag        = document.getElementById("skymag");
    mySkyMagVal     = mySkyMag.options[mySkyMag.selectedIndex].value;
    /* Query string */
    querystring += "&seeing=" + mySeeingVal + "&airmass=" + myAirmassVal + "&skymag=" + mySkyMagVal;
    

    /* Get spectro */
    mySpectro     = document.getElementById("spectro");
    mySpectroVal  = mySpectro.options[mySpectro.selectedIndex].value;
    /* myFibdiam     = document.getElementById("fibdiam");
    myFibdiamVal  = myFibdiam.options[myFibdiam.selectedIndex].value; */
    /* Query string */
    /* querystring += "&spectro=" + mySpectroVal + "&fibdiam=" + myFibdiamVal; */
    querystring += "&spectro=" + mySpectroVal + "&fibdiam=1";
    
    /* Get binning setup */
    /* mySpatBin     = document.getElementById("spatbin");
    mySpatBinVal  = mySpatBin.options[mySpatBin.selectedIndex].value;
    mySpecBin     = document.getElementById("specbin");
    mySpecBinVal  = mySpecBin.options[mySpecBin.selectedIndex].value; */
    /* Querystring */
    querystring += "&spatbin=2" + "&specbin=1";
    

    /* Get calculation parameters */
    myMethod      = document.getElementById("ETC_mode");
    myMethodVal   = myMethod.options[myMethod.selectedIndex].value;
    mySNR         = document.getElementById("snr");
    mySNRVal      = parseFloat(mySNR.value);
    myTime        = document.getElementById("etime");
    myTimeVal     = parseFloat(myTime.value);
    /* Querystring */
    querystring += "&meth=" + myMethodVal + "&etime=" + myTimeVal + "&snr=" + mySNRVal;
    
    
    /* Get target parameters */
    Type        = document.getElementById("src_type");
    TypeVal     = Type.options[Type.selectedIndex].value;
    Redshift    = document.getElementById("redshift");
    RedshiftVal = parseFloat(Redshift.value);
    TgtMag         = document.getElementById("magnitude");
    TgtMagVal      = parseFloat(TgtMag.value);
    Band      = document.getElementById("band");
    BandVal   = Band.options[Band.selectedIndex].value;
    Template    = document.getElementById("template");
    TemplateVal = Template.options[Template.selectedIndex].value;
    /* Querystring */
    querystring += "&src_type=" + TypeVal + "&tgtmag=" + TgtMagVal + "&redshift=" + RedshiftVal + "&band=" + BandVal + "&template=" + TemplateVal;

    
    /* Test values */
    if (!(isNumber(myTimeVal)) || myTimeVal<0 ) {
	alert("The exposure time for SNR calculation value is invalid, please specify a positive number.");
	error_qs++;
    }
    if (!(isNumber(mySNRVal)) || mySNRVal<0 ) {
	alert("The SNR for exposure time calculation value is invalid, please specify a positive number.");
	error_qs++;
    }
    
}


function submit()
{
    
    /* Query all parameters (with 50ms delay for safety) */
    setTimeout(updateQuerystring('submit'),50);
    
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
	
	/* Call to mse_wrapper.py */
	ajaxRequest.open("GET","../cgi-bin/mse/mse_wrapper.py" + querystring, true);
	ajaxRequest.onreadystatechange = function() {
	    if (ajaxRequest.readyState == 4) {


		/* Printing results values */
		response = ajaxRequest.responseText;
		lines = response.split("\n");
		/* Line with the parameters */
		paramtext = lines[0] + "<br>";
		/* Div tag from bokeh plot*/
		divtext = lines[2] + lines[3] + lines[4] + "<br>";
		/* Script */
		scripttext = "";
		for (i=7; i < lines.length-2; i++) {
		    scripttext = scripttext.concat(lines[i], "\n")
		}

		/* All text */
		alltext = "";
		for (i=1; i < lines.length; i++) {
		    alltext += lines[i]
		}
		
		cellParam = document.getElementById("results");
		cellParam.innerHTML = divtext;

		eval(scripttext);
		
		/* Restore appearance of Submit button */
		document.getElementById("submit_color").bgColor = "#FF0000";
		document.getElementById("submit").value = 'Submit';
		document.getElementById("submit").disabled = false;

	    } 
	}
	ajaxRequest.send(null);
    }
}
