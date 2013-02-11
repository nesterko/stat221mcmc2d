// fixing the bind() function for browsers that don't support it, to improve 
// iPad compatibility
if (!Function.prototype.bind) {
  Function.prototype.bind = function (oThis) {
    if (typeof this !== "function") {
      // closest thing possible to the ECMAScript 5 internal IsCallable function
      throw new TypeError("Function.prototype.bind - what is trying to be bound is not callable");
    }

    var aArgs = Array.prototype.slice.call(arguments, 1), 
        fToBind = this, 
        fNOP = function () {},
        fBound = function () {
          return fToBind.apply(this instanceof fNOP
                                 ? this
                                 : oThis,
                               aArgs.concat(Array.prototype.slice.call(arguments)));
        };

    fNOP.prototype = this.prototype;
    fBound.prototype = new fNOP();

    return fBound;
  };
}

function getWindowHeight() {
  var windowHeight = 0;
  if (typeof(window.innerHeight) == 'number') {
    windowHeight = window.innerHeight;
  } else {
    if (document.documentElement && document.documentElement.clientHeight) {
      windowHeight = document.documentElement.clientHeight;
    } else {
      if (document.body && document.body.clientHeight) {
	windowHeight = document.body.clientHeight;
      }
    }
  }
  return windowHeight;
}

function setFooter() {
  if (document.getElementById) {
    var windowHeight = getWindowHeight();
    if (windowHeight > 0) {
      var el = document.getElementById('headerandcontent')
        if (el) {
          var contentHeight = el.offsetHeight;
          var footerElement = document.getElementById('footer');
          var footerHeight  = footerElement.offsetHeight;
          if (windowHeight - (contentHeight + footerHeight) >= 0) {
            footerElement.style.top = ((windowHeight - (contentHeight + footerHeight)) 
                + 'px');
          } else {
            footerElement.style.top = '0px';
          }
        }
    }
  }
}

function getURIargs() {
  // Get a reference to the embedded SCRIPT tag.
  // SCRIPT tags are loaded and parsed serially.
  // So the last script loaded is the current script being parsed.  
  var all_script_tags = document.getElementsByTagName('script');
  var script_tag = all_script_tags[all_script_tags.length - 1];
  // Get the query string from the embedded SCRIPT tag's src attribute
  var query = script_tag.src.replace(/^[^\?]+\??/,''); 
  var re = new RegExp(/(^.*\/).*/)
  var path = re.exec(script_tag.src)[1];
  // Parse query string into arguments/parameter
  var vars = query.split("&");
  var args = {};
  for (var i=0;i<vars.length;i++) {
    var pair = vars[i].split("=");
    args[pair[0]] = decodeURI(pair[1]).replace(/\+/g, ' ');  // decodeURI doesn't expand "+" to a space
  }
  args.jspath = path;
  // Output detected arguments
  return args;
}

/* a helper to size the container, uses jquery */
function computeDims(containerId, w, h) {
  if (!w || !h) {
    throw("width or height arguments were not provided!");
  }
  // convert width and height to numeric
  w = +w;
  h = +h;
  // define variables for end result, absolute pixel dimensions
  var width, height;
  if (w <= 1) {
    width = w * $("#" + containerId).parent().width();
  } else {
    width = w;
  }
  if (h <= 1) {
    height = h * $(window).height();
  } else {
    height = h;
  }
  return {w : width, h : height};
}

function addScriptToHead(url, callback)
{
  var scr = document.createElement('script');
  scr.src = url;
  scr.onload = scr.onreadystatechange = function(){
    if( scr.readyState ){
      if(scr.readyState === 'complete' || scr.readyState === 'loaded'){
        scr.onreadystatechange = null;                                                  
        if (callback) callback();
      }
    } 
    else{                               
      if (callback) callback();
    }
  };
  $('head')[0].appendChild(scr);
  return 0;
}

window.onload = function() {
  setFooter();
}
window.onresize = function() {
  setFooter();
}

