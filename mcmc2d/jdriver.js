// a sample d3 visualization for Stat 221
// Metropolis Algorithm
// Konstantin Kashin, reusing some of Sergyi's code
// February 18, 2013

(function() {
  var args = getURIargs();

  // put in the custom style css for the visualization
  $('head').append('<link rel="stylesheet" href="' + 
    args.jspath + 'style.css">');
    
// Begin visualization function         
  function visual ()
{
  // this line is needed for MathJax to compile correctly
  MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
   
  // set global variables for target bivariate normal distribution
  var targetMuX = 0, targetMuY = 0, targetSigmaX=1, targetSigmaY=1, targetRho = 0.5,
  resolution = 250;

  // set global variables for proposal bivariate distribution 
  var proposalMuX = 0, proposalMuY = 0, proposalSigmaX=1, proposalSigmaY=1, proposalRho = 0,
  resolution = 250;

  // compute the total dimensions of the plot as provided by the 
  // URI arguments (computeDims lives in the *utils.js file)
  var dims = computeDims(args.c, args.w, args.h);

  var w = dims.w, h = dims.h;

// define plot margins for the plots (same for all plots)
  var m = {t : 10, r : 20, b : 40, l : 60};

  // define individual svg container width, and also the width and height
  // of the plotting g elements within each svg
  var wv = w - 10, pw = wv - m.r - m.l, ph = h - m.t - m.b;

  function init (container, w, h, title, xlab, ylab, data)
  {
    // initialize the SVG container (this will contain the contour plot)
    var con = d3.select("#" + container).insert("div", ":first-child")
      .attr('class', 'ivcont');

    // add visualization title
    con.append("h2").html(title);

    // initialize the SVG element
    var svg = con.append("svg")
      .attr("width", w)
      .attr("height", h);

    // add a clip path so that nothing outside the clip path gets plotted
    var clip = svg.append("defs").append("clipPath")
      .attr("id", "clip")
      .append("rect")
      .attr("width", pw)
      .attr("height", ph);

    // initialize the axes components
    var x = d3.scale.linear().range([0, pw]),
        y = d3.scale.linear().range([ph, 0]),
        xAxis = d3.svg.axis().scale(x).tickSize(3, 3, 1),
        yAxis = d3.svg.axis().scale(y).tickSize(3, 3, 1).orient("left");

    // initialize the domains of x and y
    // the domain of x will be the same throughout : [0, 1], but y's will change
    x.domain([0, 1]);
    y.domain([0, 1]);

    // initialize the g elements inside the svg container for x and y 
    // axes, and also for the plot itself. We are also calling the 
    // axis functions on the containers to keep track of domain changes, 
    // d3 shines here

    var ycon = svg.append("g")
      .attr("class", "y axis")
      .attr("transform", "translate(" + m.l + "," + m.t + ")")
      .call(yAxis);

    var xcon = svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(" + m.l + "," + parseFloat(ph + m.t) + ")")
      .call(xAxis);

    // add the axes labels (if there are any)
    if (xlab) {
      xcon.append("text").attr("class", "label").attr("x", pw / 2)
        .attr("y", 0).attr("dy", "+2em").attr("text-anchor", "middle")
        .text(xlab);
    }
    if (ylab) {
      ycon.append("text").attr("class", "label").attr("x", 0).attr("y", 0)
        .attr("dy", "-2.5em").attr("text-anchor", "middle")
        .text(ylab)
        .attr("transform", "rotate(-90) translate(-" +
            parseFloat(ph / 2) + "," + 0 + ")");
    }

    // the actual plotting area element within the plots
    var plot = svg.append("g")
      .attr("transform", "translate(" + m.l + "," + m.t + ")");

    // return a handy Object with all the needed elements
    return {con : con, svg : svg, plot : plot, data : data,
      ycon : ycon, xcon : xcon,
      x : x, y : y, xAxis : xAxis, yAxis : yAxis};
  }

  function redraw ()
  {
    // draw contour plot (fxn defined below)
    drawContour(contour);
    return 0;
  }

  // generate data for the plots
  function generate ()
  {
   // create data for bivariate normal contour plot
  	data.contour = makeBivarNormContourData(targetMuX,targetMuY,targetSigmaX,targetSigmaY,targetRho,resolution,nLevels,cThreshold)
    
    // initialize Metropolis algorithm
    metropolisInit();

    return 0;
  }



  /*************************** HELPER FXNS *********************************/

  // Contour plot settings
  var nLevels = 12, cThreshold = 0.05;

  // Define density of bivariate normal distribution
  function bivariateNormal(x,y,muX,muY,sigmaX,sigmaY,rho){
 		var det = Math.pow(sigmaX,2)*Math.pow(sigmaY,2)-Math.pow(rho,2)*Math.pow(sigmaX,2)*Math.pow(sigmaY,2);
 		var epower = -1/2*1/(1-Math.pow(rho,2))*(Math.pow((x-muX),2)/Math.pow(sigmaX,2) - 2*rho*(x-muX)*(y-muY)/(sigmaX*sigmaY)+Math.pow((y-muY),2)/Math.pow(sigmaY,2));
 		return 1/(2*Math.PI)*Math.pow(det,-1/2)*Math.exp(epower);
 	}
 
  // This is a function to make a contour plot of a bivariate normal
  // Note that the actual densities are calculated within the fxn
  function makeBivarNormContourData (muX,muY,sigmaX,sigmaY,rho,resolution,nlev, thresh)
   {  	
  	// calculate ranges 
  	var xRange = [muX-3*sigmaX, muX+3*sigmaX];
  	var yRange = [muY-3*sigmaY, muY+3*sigmaY];
  	    
    // create xx and yy arrays (which give support of density)
    var xx = jStat.seq(xRange[0], xRange[1], resolution);
    var yy = jStat.seq(yRange[0], yRange[1], resolution);
   
   	// create nested array of densities (first level is values of xx, second level is values of yy)
   	var sur = jStat.seq(0, xx.length - 1, xx.length).map(function (jj) 
        {
          return jStat.seq(0, yy.length - 1, yy.length).map(function (ii)
            {
                return bivariateNormal(xx[jj],yy[ii],muX,muY,sigmaX,sigmaY,rho);
            });
        });
    
    // tease out all available z values
    var zz = [];
    $M(sur).map(function (e) { zz.push(e); return 0;});
    // determine the levels
    var ra = range(zz), ad = (ra[1] - ra[0]) / (nlev - 1);
    var levels = jStat.seq(0, nlev - 1, nlev).map(function (ii)
        {return ra[0] + ad * ii; });

    // use Jason Davies' implementation of the contour plot
    var c = new Conrec();
    var co = c.contour(sur, 0, xx.length - 1, 0, yy.length - 1,
        xx, yy, levels.length, levels);

    return {xx : xx, yy : yy, points : c.contourList()};
  }
  

   // This is a function to make a contour plot from two univariate distributions
   function makeContourData (xdat, ydat, nlev, thresh)
  {
    // calculate the values at the grid
    var xx = xdat.xx, yy = ydat.xx;
    var sur = jStat.seq(0, xx.length - 1, xx.length).map(function (jj) 
        {
          return jStat.seq(0, yy.length - 1, yy.length).map(function (ii)
            {
                return xdat.yy[jj] * ydat.yy[ii];
            });
        });
    // tease out all available z values
    var zz = [];
    $M(sur).map(function (e) { zz.push(e); return 0;});
    // determine the levels
    var ra = range(zz), ad = (ra[1] - ra[0]) / (nlev - 1);
    var levels = jStat.seq(0, nlev - 1, nlev).map(function (ii)
        {return ra[0] + ad * ii; });

    // use Jason Davies' implementation of the contour plot
    var c = new Conrec();
    var co = c.contour(sur, 0, xx.length - 1, 0, yy.length - 1,
        xx, yy, levels.length, levels);

    return {xx : xdat.xx, yy : ydat.xx, points : c.contourList()};
  }
  

  // Function to create array of pair objects (some point in support of density and then the density at that point)
  function makeDensityData (ran, resolution, densityFn, normalize)
  {
    var tmp = Object();
    tmp.xx = jStat.seq(ran[0], ran[1], resolution); // create array of 250 values within range of RV
    // store associated densities in yy array
    tmp.yy = tmp.xx.map(function (e)
        { return densityFn(e); });
    // max and min densities
    var maxdensity = d3.max(tmp.yy);
    var mindensity = d3.min(tmp.yy);
    // rescale make maxdensity = 1
    if (normalize) {
      tmp.yy = tmp.yy.map(function (e) { return e / maxdensity; });
    }
    tmp.maxdensity = maxdensity;
    tmp.mindensity = mindensity;
    // make array of pair objects
    tmp.points = makePoints(tmp.xx, tmp.yy);
    return tmp;
  }
  function normalizeDensities ()
  {
    // get the overall maximum value of the density
    var maxx = d3.max(data.den.map(function (el)
          {return el.data.maxdensity;}));
    data.maxden = maxx;
    // normalize density values by dividing each value by the overall maximum
    for (var ii = 0; ii < data.den.length; ii++) {
      var te = data.den[ii].data.yy.map(function (el) {return el / maxx;});
      data.den[ii].data.maxdensity = d3.max(te);
      data.den[ii].data.mindensity = d3.min(te);
      data.den[ii].data.yy = te;
      data.den[ii].data.points = makePoints(data.den[ii].data.xx, te);
    }
    return 0;
  }
  // function to generate normal density
  function normal (x, mu, ss)
  {
    var fun = jStat.normal.pdf;
    return fun(x, mu, ss);
  }

  // function to create array of pair objects (points) from two arrays
  function makePoints (xx, yy)
  {
    return jStat.seq(0, xx.length - 1, xx.length).map (function (ii)
        {
          return { x : xx[ii], y : yy[ii] };
        });
  }
  function range (ar) { return [d3.min(ar), d3.max(ar)]; }
  
  /*************************** METROPOLIS ALGORITHM *********************************/
  
  // function to initialize Metropolis algorithm
  function metropolisInit ()
  {
    data.metropolis = Object();
    data.metropolis.x = Object();
    data.metropolis.y = Object();
	// define initial values
    data.metropolis.x.val = targetMuX;
    data.metropolis.y.val = targetMuY;
    // define arrays that will store all values
    data.metropolis.x.vals = [];
    data.metropolis.y.vals = [];
    // define array that will store if proposal was accepted at each iteration
    data.metropolis.accepted = [];
    data.metropolis.currentIteration = 0;

	// set acceptance rate in legend equal to blank (useful for when "refresh" plot)    
    var acceptRate = d3.selectAll(".acceptRate");
	acceptRate.html("Acceptance Rate of Proposals: - %");;
		
    return 0;
  }
  
 // define function to do next iteration
 function metropolisNextIteration(){
  	metropolisJump();
 	redraw();
 }
 // define function to do next 100 iterations
 function metropolis100Iterations(){
  	for (var ii=0; ii<100; ii++) metropolisJump();
 	redraw();
 }
 // define function to do next 1000 iterations
 function metropolis1000Iterations(){
 	for (var ii=0; ii<1000; ii++) metropolisJump();
 	redraw();
 } 
 
 // define function to reset sampling
 function metropolisReset () {metropolisInit(); redraw;}
 
 // define function to do acceptance/rejection jump
 function metropolisJump()
 {
 	// define current position (this is X_t,Y_t)
 	var currentX  = data.metropolis.x.val; 
 	var currentY =  data.metropolis.y.val;
 	
 	// proposal step: sample candidate X,Y from proposal distribution
 	// proposal distribution here is bivariate normal with variance and rho parameters
 	// defined by the user
 	var xRaw = jStat.normal.sample(0,1);
 	var xProposal = proposalSigmaX*xRaw + currentX; 
 	var yRaw = jStat.normal.sample(0,1);
 	var yProposal = proposalSigmaY*(proposalRho*xRaw +Math.sqrt(1-Math.pow(proposalRho,2))*yRaw) + currentY
 	
 	
 	// calculate probability of acceptance as ratio of densities evaluated at proposal to original point
 	var densOrig = bivariateNormal(currentX,currentY,targetMuX,targetMuY,targetSigmaX,targetSigmaY,proposalRho);
 	var densProposal = bivariateNormal(xProposal,yProposal,targetMuX,targetMuY,targetSigmaX,targetSigmaY,proposalRho);
 	
 	var probAccept = densProposal/densOrig;
 	
 	// bound acceptance probability at a max of 1
 	var probAcceptBounded = Math.min(densProposal,densOrig)
 
 	// select next step (accept or reject proposal)
 	unifIter = jStat.uniform.sample(0,1);
 	// Accept with probability of probAcceptBounded
 	if(unifIter < probAcceptBounded){
 		data.metropolis.accepted.push(true);
 		data.metropolis.x.val = xProposal;
 		data.metropolis.y.val = yProposal;
 		data.metropolis.x.vals.push(xProposal);
 		data.metropolis.y.vals.push(yProposal); 
 	} else{
 	// reject proposal
 	data.metropolis.accepted.push(false);
 	data.metropolis.x.val = currentX;
 	data.metropolis.y.val = currentY;
 	data.metropolis.x.vals.push(currentX);
 	data.metropolis.y.vals.push(currentY); 
 	}
 	data.metropolis.currentIteration++;
  } // end of metropolisJump function

 
 
 

/******************** INITIALIZE DATA OBJECT ******************/
  // generate/initialize the data object
  var data = Object();
  generate();

  // initialize data object & all data for plots
  var contour = init(args.c,wv,h, "Metropolis Sampling from Bivariate Normal", "X", "Y", data)

  
/********************** DEFINE LEGEND ***********************/


  // a function to add the legend
  function legend ()
  {
    // select where we'll be adding the legend to
    var leg = d3.select("#" + args.c).append("div")
     .attr('class', "ivcont userviscontrols");
    // set the width of the container
    // add all fields, text inputs and the button
    
    // Create buttons to control iterations of Metropolis algorithm
 	// Define what name of button is and the fxn it calls when pressed
    var dataIterBut = [
      {name : "Next Iteration", fn : metropolisNextIteration},
      {name : "100 Iterations", fn : metropolis100Iterations},
      {name : "1000 Iterations", fn : metropolis1000Iterations}
      //{name : "Refresh", fn : metropolisReset}
    ];
    var iterControls = leg.append('div').attr('class','iterControls');
    // Create buttons
    var iterButtons = iterControls
      .selectAll(".button").data(dataIterBut);
    iterButtons.enter().append("input")
      .attr("type", "submit").attr("class", "button")
      .attr("value", function (d) { return d.name; })
      .on("click", function (d) { return d.fn(); });
      
    // Create field to display acceptance rate of proposals in legend
    var acceptRate = iterControls.append('h3').attr('class','acceptRate');
	acceptRate.html("Acceptance Rate of Proposals: - %")
      
    // Create fields to input parameters for target and proposal distributions
    // Define input data for target distribution (this will be different input blanks) 
    var inputdataTarget = [
    {lab : '\\\\( \\mu_X \\\\)', cls : 'targetMuX', val : targetMuX},
    {lab : '\\\\( \\mu_Y \\\\)', cls : 'targetMuY', val : targetMuY},
    {lab : '\\\\( \\sigma_X \\\\)', cls : 'targetSigmaX', val : targetSigmaX},
    {lab : '\\\\( \\sigma_Y \\\\)', cls : 'targetSigmaY', val : targetSigmaY},
    {lab : '\\\\( \\rho \\\\)', cls : 'targetRho', val : targetRho}];
    
    // Define input data for proposal distribution (this will be different input blanks) 
    var inputdataProposal = [
    //{lab : '\\\\( \\mu_X \\\\)', cls : 'proposalMuX', val : proposalMuX},
    //{lab : '\\\\( \\mu_Y \\\\)', cls : 'proposalMuY', val : proposalMuY},
    {lab : '\\\\( \\sigma_X \\\\)', cls : 'proposalSigmaX', val : proposalSigmaX},
    {lab : '\\\\( \\sigma_Y \\\\)', cls : 'proposalSigmaY', val : proposalSigmaY},
    {lab : '\\\\( \\rho \\\\)', cls : 'proposalRho', val : proposalRho}];
    
    // Append div for target distribution and populate
    var tdist = leg.append('div').attr('class','targetdist');
    tdist.append('h3').html('Target Distribution');
    var inpTarget = tdist.selectAll(".labsandpars").data(inputdataTarget)
      .enter().append('div').attr('class', 'labsandpars');
    // this creates the labels for the input boxes
    inpTarget.insert('p').html(function (d) {return d.lab; });
    // This create the input boxes
    inpTarget.append('input').attr('class', function (d) { return 'pars ' + d.cls})
      .attr('type', 'number')
      .attr('value', function (d) {return d.val;})
      .attr('step', function (d) {
      if(d.cls=="targetRho"){
      return 0.1; // steps of 0.1 for rho parameter
      }
      else{return 1;} // otherwise steps of 1
      })
      .attr('min', function (d) {
      if(d.cls=="targetRho"|d.cls=="targetSigmaX"|d.cls=="targetSigmaY"){
      return 0; // minimum value of 0 for rho and SDs
      }
      })
      .attr('max', function (d) {
      if(d.cls=="targetRho"){
      return 1; // max value of 1 for rho
      }
      });
    
    // Append div for proposal distribution and populate
    var pdist = leg.append('div').attr('class','propdist');
    pdist.append('h3').html('Proposal Distribution');
    var inpProposal = pdist.selectAll(".labsandpars").data(inputdataProposal)
      .enter().append('div').attr('class', 'labsandpars');
    // this creates the labels for the input boxes
    inpProposal.insert('p').html(function (d) {return d.lab; });
    // this creates the input boxes
    inpProposal.append('input').attr('class', function (d) { return 'pars ' + d.cls})
      .attr('type', 'number')
      .attr('format', function (d, i) 
          { return i == 1? '\\d+' : '\\d'; })
      .attr('value', function (d) {return d.val;})
      .attr('step', function (d) {
      if(d.cls=="proposalRho"){
      return 0.1;
      }
      else{return 1;}
      })
      .attr('min', function (d) {
      if(d.cls=="proposalRho"|d.cls=="proposalSigmaX"|d.cls=="proposalSigmaY"){
      return 0;
      }
      })
      .attr('max', function (d) {
      if(d.cls=="proposalRho"){
      return 1;
      }
      });
      
    // Create button to refresh / reset plot
    // Define what name of button is and the fxn it calls when pressed
    var databut = [
      {name : "Refresh", fn : getPars}
    ];
    
    // This creates refresh button
    var buttons = leg.append('div').attr('class','refresh')
        .selectAll(".refreshButton").data(databut);
    buttons.enter().append("input")
      .attr("type", "submit").attr("class", "refreshButton")
      .attr("value", function (d) { return d.name; })
      .on("click", function (d) { return d.fn(); });
    return 0;
  }

  // Define fxn that runs when hit "Refresh"
  function getPars ()
  {
  // extract values of all input forms of class "pars"
  var boxes = d3.selectAll("#" + args.c + " .pars");
  var vals = boxes[0].map(
      function (el) {return +d3.select(el).property("value");});
  // update global values of parameters with new user-specified values
  targetMuX = vals[0]; targetMuY = vals[1]; targetSigmaX = Math.abs(vals[2]);
  targetSigmaY = Math.abs(vals[3]); targetRho = vals[4];
  // Note: need to do same here for proposal distribution (positions 5-9 of vals)!!!
	proposalSigmaX = vals[5]; proposalSigmaY = vals[6]; proposalRho = vals[7];
 
  // regenerate data
  generate();
  // redraw svg
  redraw();
  return 0;
  }

 
  function bound (minmax, ax)
  {
    var mm = minmax;
    var su = data.den;
    var va = d3[mm](su.map(function (el) {return d3[mm](el.data[ax]);}));
    return va;
  }

  function getBounds (ax) 
  {
    return ['min', 'max'].map( function (mm)
        { return bound(mm, ax); });
  }



  /************************** DRAWING FXNS *********************************/
  
  
  // This function draws the contour plot + axes
  function drawContour (obj)
  {
    // set the axis right
    obj.y.domain(range(obj.data.contour.yy));
    obj.x.domain(range(obj.data.contour.xx));

    // draw the target density
    var line = d3.svg.line()
      .x(function(d, i) { return obj.x(d.x); })
      .y(function(d, i) { return obj.y(d.y); });

	// obj.data.contour.points is array where element is an array corresponding to each contour line
    var lines = obj.plot.selectAll("path.contour")
      .data(obj.data.contour.points);
     
    // This draws the path for the contour plot
    // Note that input is line, which array with each element corresponding to a curve on the contour plot
    lines.enter().append("path")
      .attr("class", "contour")
      .attr("d", line)

    // draw the point at which the Metropolis algorithm is at
    var poi = {x : obj.data.metropolis.x.val, y : obj.data.metropolis.y.val};

    var point = obj.plot.selectAll("circle")
      .data([poi]);

    point.enter().append("circle")
      .attr("class", "currentValue")
      .attr("cx", function (d, i) { return obj.x(d.x); })
      .attr("cy", function (d, i) { return obj.y(d.y); })
      .attr("r", 7)
      .each( function (d)
          { this._current = Object(); this._current.d = d;
            this._current.name = 'contour'; this._current.obj = obj; return 0; });


	// draw previous points (the entire array of sampled points, including current point)     
	// create pairs of points
	 var prevPoi = makePoints(data.metropolis.x.vals,data.metropolis.y.vals)
	// create circles of class "prev" and pass data
    var prevPoints = obj.plot.selectAll("circle.prev")
      .data(prevPoi);
	// plot points
    prevPoints.enter().append("circle")
      .attr("class", "prev")
      .attr("cx", function (d, i) { return obj.x(d.x); })
      .attr("cy", function (d, i) { return obj.y(d.y); })
      .attr("r", 4)


    // create object that will hold lines connecting sampled points; plot lines
    var connectLines = []
    if(data.metropolis.currentIteration > 0){
   		for(var ii=1; ii < prevPoi.length; ii++){
    	connectLines.push({'x1':prevPoi[ii-1].x,'y1':prevPoi[ii-1].y,'x2':prevPoi[ii].x,'y2':prevPoi[ii].y})
    	}
    }

    var lines = obj.plot.selectAll("line")
    .data(connectLines);
    
  lines.enter().append("line")
    .attr('class', 'connect')
    .attr("x1", function (d, i) { return obj.x(d.x1); })
    .attr("y1", function (d, i) { return obj.y(d.y1); })
     .attr("x2", function (d, i) { return obj.x(d.x2); })
    .attr("y2", function (d, i) { return obj.y(d.y2); })
    .style("opacity", 0.85)
    .style("stroke", "red")
    .attr("stroke-dasharray", "5, 5");
  
    

	// Calculate acceptance rate of Metropolis algorithm (only if at least 1 iteration...)
	if(data.metropolis.accepted.length!=0)
	{
		var sumAccept = data.metropolis.accepted.reduce(function(a, b) { return a + b });
		var percAccept = sumAccept/data.metropolis.accepted.length;
		var acceptRate = d3.selectAll(".acceptRate");
		acceptRate.html("Acceptance Rate of Proposals: " + percAccept.toFixed(5)+"%");;
	}


    // handle the transitions (this is necessary to reset plot when hit "Refresh")
    var tr = obj.svg.transition().duration(500);
    // rescale axes
    tr.select(".y.axis").call(obj.yAxis);
    tr.select(".x.axis").call(obj.xAxis);
    // Redraw current point in Metropolis algorithm + path leading to it
    tr.select("circle").attr("cx", function(d){return obj.x(d.x);})
      .attr("cy",function(d){return obj.y(d.y);});
    // Redraw contour plot
    tr.selectAll("path.contour").attr('d', function (d)
       { return line(d) } );
    
    // transitions for lines
    tr.selectAll("line.connect")
    .attr("x1", function (d, i) { return obj.x(d.x1); })
    .attr("y1", function (d, i) { return obj.y(d.y1); })
     .attr("x2", function (d, i) { return obj.x(d.x2); })
    .attr("y2", function (d, i) { return obj.y(d.y2); })
    .style("opacity", 0.85)
    .style("stroke", "red")
    .attr("stroke-dasharray", "5, 5");
    
    // remove all previous points and lines connecting them from plot when hit "refresh"
    prevPoints.exit().remove();
    lines.exit().remove();
    return 0;
  }

 /*************************** INITIALIZE PLOT ******************************/
  legend();
  redraw();
  getPars();
  
  
  // register an event listener for pressing the Enter key when editing numbers
  d3.selectAll("#" + args.c + " .pars").on("keyup", function ()
      {
      if (d3.event.keyCode == 13) getPars();
    });

  return 0;
} // end of visual() fxn

// call visual function
visual();
//addScriptToHead(args.jspath+"conrec.js",visual);

return 0;
})();
