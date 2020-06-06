/**
 * Attempt at using AJAX to "submit form and then plot SNR spectrum"
 * @author Nicolas FLAGEY, Billy Mahoney
 */


function initFullFormAjaxUpload() {
    var form = document.getElementById('form_id');
    form.onsubmit = function() {
	// FormData receives the whole form
	var formData = new FormData(form);

	// We send the data where the form wanted
	var action = '/cgi-bin/mse/mse_alloc_wrapper.py';

	// Code common to both variants
	sendXHRequest(formData, action);

	// Avoid normal form submission
	return false;
    }
}


// Once the FormData instance is ready and we know
// where to send the data, the code is the same
// for both variants of this technique
function sendXHRequest(formData, uri) {
    
    /* Change appearance of Submit button */
    document.getElementById("submit").disabled = false;
    document.getElementById("submit_color").bgColor = "#FF0000";
    document.getElementById("submit").value = 'Processing request ...';
    document.getElementById("submit").disabled = true;
    
    
    // Get an XMLHttpRequest instance
    var xhr = new XMLHttpRequest();

    // Set up events
    xhr.upload.addEventListener('loadstart', onloadstartHandler, false);
    xhr.upload.addEventListener('progress', onprogressHandler, false);
    xhr.upload.addEventListener('load', onloadHandler, false);
    xhr.addEventListener('readystatechange', onreadystatechangeHandler, false);

    // Set up request
    xhr.open('POST', uri, true);

    // Fire!
    xhr.send(formData);

    /* Restore appearance of Submit button */
    /* document.getElementById("submit_color").bgColor = "#FF0000";
    document.getElementById("submit").value = 'Submit';
    document.getElementById("submit").disabled = false; */
    
}

// Handle the start of the transmission
function onloadstartHandler(evt) {
    var div = document.getElementById('upload-status');
    div.innerHTML = 'Upload started.';
}

// Handle the end of the transmission
function onloadHandler(evt) {
    var div = document.getElementById('upload-status');
    div.innerHTML += '<' + 'br>File uploaded. Waiting for response.';
}

// Handle the progress
function onprogressHandler(evt) {
    var div = document.getElementById('progress');
    var percent = evt.loaded/evt.total*100;
    div.innerHTML = 'Progress: ' + percent + '%';
}

// Handle the response from the server
function onreadystatechangeHandler(evt) {
    var status, text, readyState;

    try {
	readyState = evt.target.readyState;
	text = evt.target.responseText;
	status = evt.target.status;
    }
    catch(e) {
	return;
    }

    if (readyState == 4 && status == '200' && evt.target.responseText) {
	var status = document.getElementById('upload-status');
	var result = document.getElementById('results');

	if (evt.target.responseText == 'Your file does not have the right format.') {
	    status.innerHTML += 'Your file does not have the right format.';
	    result.innerHTML = 'Please check the format of your file and re-upload it before submitting.';
	} else {
	    status.innerHTML += '<' + 'br>Success!';
	    result.innerHTML = evt.target.responseText;	    
	}
    }
}
