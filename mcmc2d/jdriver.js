// a sample d3 visualization for Stat 221
// Strength of prior
// sergiy 20120208


(function() {
  var args = getURIargs();

  // put in the custom style css for the visualization
  $('head').append('<link rel="stylesheet" href="' + 
    args.jspath + 'style.css">');
        
  function visual ()
{
  // global variables for likelihood and prior means and variances
 // var xbar = 0, ss = 2, mu0 = 1, ss0 = 1;
  // a global variable to tell the number of observations
 // var nn = 10;

  // compute the total dimensions of the plot as provided by the 
  // URI arguments (computeDims lives in the *utils.js file)
  var dims = computeDims(args.c, args.w, args.h);

  var w = dims.w,
h = dims.h;

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
    //drawDen(den);
    drawAug(aug);
    return 0;
  }

  // generate data for the plots
  function generate ()
  {
  /*  getRanges();
    data.den = ['prior', 'likelihood', 'posterior'].map(function (nam)
        {
          var ran, res, fn;
          if (nam == 'prior') {
            ran = priRange; res = priResolution; fn = prior;
          } else if (nam == 'likelihood') {
            ran = likRange; res = likResolution; fn = likelihood;
          } else if (nam == 'posterior') {
            ran = postRange; res = postResolution; fn = posterior;
          }
          var dat = makeDensityData(ran, res, fn, false);
          return { name : nam, data : dat};
        });
    //normalizeDensities();
    window.data = data;
    */
    // calculate target density data
    data.x = makeDensityData(xRange, xResolution, xDensity, true);
    data.y = makeDensityData(yRange, yResolution, yDensity, true);
    //data.aug = makeContourData(data.x, data.y, nLevels, cThreshold);
	data.aug = makeBivarNormContourData(targetMuX,targetMuY,targetSigmaX,targetSigmaY,targetRho,resolution,nLevels,cThreshold)
    // initialize Metropolis
    metropolisInit();

    return 0;
  }

 /* // variables for resolutions and ranges of the likelihood, prior, and 
  // posterios
  var likResolution = 250, priResolution = 250, postResolution = 250,
      likRange, priRange, postRange;

  function getRanges ()
  {
    var nsd = 3.5;
    var lsd = ss / Math.sqrt(nn);
    likRange = [xbar - nsd * lsd, xbar + nsd * lsd];
    priRange = [mu0 - nsd * ss0, mu0 + nsd * ss0];
    var pp = posteriorMeanSD();
    postRange = [pp.mean - nsd * pp.sd, pp.mean + nsd * pp.sd];
    return 0;
  }
  */
  
  /***** HELPER FXNS *********************************/

  // nuts and bolts functions - simulations, HMC steps, contour plots etc
  // numerics related constants
  //var tobs = [-5, 5];
  //var posRange = [-9, 9], posResolution = 250;
  //var cauchyscale = 0.5;
  //var momMean = 0, momSD = 1, momResolution = 250,
  //    momRange = [momMean - 3 * momSD, momMean + 3 * momSD];
  // Contour plot settings
  var nLevels = 12, cThreshold = 0.05;

  var targetMuX = 0;
  var targetMuY = 0;
  var targetSigmaX=1;
  var targetSigmaY=1;
  var targetRho = 0.5;
  var resolution = 250;
  
  var proposalMuX = 0;
  var proposalMuY = 0;
  var proposalSigmaX=1;
  var proposalSigmaY=1;
  var proposalRho = 0;
  var resolution = 250;
  
 
  // in the new simulation, we will have bivariate normal as target distribution
  var xMean = 0;
  var xSD = 2;
  var xRange = [xMean-3*xSD, xMean+3*xSD]; // calculated as -6,6
  var xResolution=250;
  
  var yMean = 0;
  var ySD = 1;
  var yRange = [yMean-3*ySD, yMean+3*ySD]; // calculated as -3,3
  var yResolution = 250;
  
  
  // replace these next 2 fxns with normal density below
  function xDensity (x){
  var pdf = jStat.normal.pdf;
  return pdf(x, xMean, xSD);
  }

  function yDensity (y){
  var pdf = jStat.normal.pdf;
  return pdf(y, yMean, ySD);
  }
  
  
  // what if y was exponential?
  xRate = 0.5;
  xRangeExpo = [0,5];
  
  function xDensityExpo (x){
  var pdf = jStat.exponential.pdf;
  return pdf(x, xRate);
  }
 
  
  // what if y was beta?
  xAlpha = 2;
  xBeta = 2;
  xRangeBeta = [0,1];
  
  function xDensityBeta (x){
  var pdf = jStat.beta.pdf;
  return pdf(x, xAlpha, xBeta);
  }
  
  
  // what if X was Cauchy?
  xCauchyScale = 2;
  xCauchyLocation = 0;
  xRangeCauchy = [-8,8];
  
  function xDensityCauchy (x){
  var pdf = jStat.cauchy.pdf;
  return pdf(x,xCauchyLocation, xCauchyScale);
  }
  
  // what if we had a bivariate normal?
  var xMean = 0;
  var xSD = 1;
  var xRange = [xMean-3*xSD, xMean+3*xSD]; // calculated as -6,6
  var xResolution=250;
  
  var yMean = 0;
  var ySD = 1;
  var yRange = [yMean-3*ySD, yMean+3*ySD]; // calculated as -3,3
  var yResolution = 250;
  
  var rho = 0.5;
 
 
 
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
  
 
 	function bivariateNormal(x,y,muX,muY,sigmaX,sigmaY,rho){
 		var det = Math.pow(sigmaX,2)*Math.pow(sigmaY,2)-Math.pow(rho,2)*Math.pow(sigmaX,2)*Math.pow(sigmaY,2);
 		var epower = -1/2*1/(1-Math.pow(rho,2))*(Math.pow((x-muX),2)/Math.pow(sigmaX,2) - 2*rho*(x-muX)*(y-muY)/(sigmaX*sigmaY)+Math.pow((y-muY),2)/Math.pow(sigmaY,2));
 		return 1/(2*Math.PI)*Math.pow(det,-1/2)*Math.exp(epower);
 	}
 

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
  
 
function makeBivarDensityData (ran, resolution, densityFn, normalize)
  {
  	// initialize empty object
    var tmp = Object();
    // create xx and yy arrays (which give support of density)
    tmp.xx = jStat.seq(ran[0], ran[1], resolution); // create array of 250 values within range of RV
    tmp.yy = jStat.seq(ran[0], ran[1], resolution); // create array of 250 values within range of RV

	// how to pass on two commands here!
    tmp.zz = tmp.xx.map(function (e)
        { return densityFn(e); });
        
    // max and min density
    var maxdensity = d3.max(tmp.zz);
    var mindensity = d3.min(tmp.zz);
    // normalize so that max density = 1
    if (normalize) {
      tmp.yy = tmp.yy.map(function (e) { return e / maxdensity; });
    }
    tmp.maxdensity = maxdensity;
    tmp.mindensity = mindensity;
    tmp.points = makePoints(tmp.xx, tmp.yy);
    return tmp;
  }
  
  	
 	
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
    // make maxdensity = 1
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
  
  //******** NEW **********//
  // Metropolis portion of the code
  var nsteps = 7, stepsize = 0.3;
  function metropolisInit ()
  {
    data.metropolis = Object();
    data.metropolis.x = Object();
    data.metropolis.y = Object();
    data.metropolis.x.val = data.metropolis.y.val = 0;
    data.metropolis.steps = [ {x : 0, y : 0} ];
    data.metropolis.x.vals = [];
    data.metropolis.y.vals = [];
    data.metropolis.currentStep = -1;
    data.metropolis.currentIteration = 0;
    return 0;
  }
  
  //******** END NEW **********//
  function posteriorMeanSD ()
  {
    var tt = ss * ss / nn, sq = ss0 * ss0;
    return { mean : (tt * mu0 + sq * xbar) / (tt + sq), 
      sd : Math.sqrt(Math.pow(Math.pow(tt, -1) + Math.pow(sq, -1), -1))
    }
  }
  function likelihood (par) { return normal(par, xbar, ss / Math.sqrt(nn)); }
  function prior (par) { return normal(par, mu0, ss0); }
  function posterior (par) 
  {
    var pp = posteriorMeanSD();
    return normal(par, pp.mean, pp.sd);
  }
  // function to create array of pair objects (points)
  function makePoints (xx, yy)
  {
    return jStat.seq(0, xx.length - 1, xx.length).map (function (ii)
        {
          return { x : xx[ii], y : yy[ii] };
        });
  }
  function range (ar) { return [d3.min(ar), d3.max(ar)]; }



/********************INITIALIZE PLOT******************/
  // generate/initialize the data object
  var data = Object();
  generate();

  // initialize our plots and data object
  //var den = init(args.c, wv, h, "Functions", 'parameter', 
  //    "function value", data);  
  // new
  var aug = init(args.c,wv,h, "Bivariate Normal", "X", "Y", data)

  
/******************LEGEND *********************/


  // a function to add the legend
  function legend ()
  {
    // select where we'll be adding the legend to
    var leg = d3.select("#" + args.c).append("div")
     .attr('class', "ivcont userviscontrols");
    // set the width of the container
    // add all fields, text inputs and the button
    var data = [
    {cls : 'prior', name : 'Prior'},
    {cls : 'likelihood', name : 'Likelihood'},
    {cls : 'posterior', name : 'Posterior'}
    ],
      databut = [
      {name : "Refresh", fn : getPars}
    ];
    var inputdataTarget = [
    {lab : '\\\\( \\mu_X \\\\)', cls : 'targetMuX', val : targetMuX},
    {lab : '\\\\( \\mu_Y \\\\)', cls : 'targetMuY', val : targetMuY},
    {lab : '\\\\( \\sigma_X \\\\)', cls : 'targetSigmaX', val : targetSigmaX},
    {lab : '\\\\( \\sigma_Y \\\\)', cls : 'targetSigmaY', val : targetSigmaY},
    {lab : '\\\\( \\rho \\\\)', cls : 'targetRho', val : targetRho}];
    var inputdataProposal = [
    {lab : '\\\\( \\mu_X \\\\)', cls : 'proposalMuX', val : proposalMuX},
    {lab : '\\\\( \\mu_Y \\\\)', cls : 'proposalMuY', val : proposalMuY},
    {lab : '\\\\( \\sigma_X \\\\)', cls : 'proposalSigmaX', val : proposalSigmaX},
    {lab : '\\\\( \\sigma_Y \\\\)', cls : 'proposalSigmaY', val : proposalSigmaY},
    {lab : '\\\\( \\rho \\\\)', cls : 'proposalRho', val : proposalRho}];
     var tdist = leg.append('div').attr('class','targetdist');
    tdist.append('h3').html('Target Distribution');
    var inpTarget = tdist.selectAll(".labsandpars").data(inputdataTarget)
      .enter().append('div').attr('class', 'labsandpars');
    inpTarget.insert('p').html(function (d) {return d.lab; });
    inpTarget.append('input').attr('class', function (d) { return 'pars ' + d.cls})
      .attr('type', 'number')
      .attr('value', function (d) {return d.val;})
      .attr('step', function (d) {
      if(d.cls=="targetRho"){
      return 0.1;
      }
      else{return 1;}
      })
      .attr('min', function (d) {
      if(d.cls=="targetRho"|d.cls=="targetSigmaX"|d.cls=="targetSigmaY"){
      return 0;
      }
      })
      .attr('max', function (d) {
      if(d.cls=="targetRho"){
      return 1;
      }
      });
    
    var pdist = leg.append('div').attr('class','propdist');
    pdist.append('h3').html('Proposal Distribution');
    var inpProposal = pdist.selectAll(".labsandpars").data(inputdataProposal)
      .enter().append('div').attr('class', 'labsandpars');
    inpProposal.insert('p').html(function (d) {return d.lab; });
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
    
    
    //var buttonsandmore = leg.append('div').attr('class', 'buttons');
    var buttons = leg
      .selectAll(".button").data(databut);
    buttons.enter().append("input")
      .attr("type", "submit").attr("class", "button")
      .attr("value", function (d) { return d.name; })
      .on("click", function (d) { return d.fn(); });
    /*var table = leg.append('div').attr('class', 'tabcont')
      .append("div").attr("class", "table")
      .append("div").attr("class", "trow");
    var tcell = table.selectAll(".tcell")
      .data(data);
    var cells = tcell.enter()
      .append('div').attr('class', 'tcell');
    cells.append("svg")
      .style("display", "inline-block")
      .attr("width", 50).attr("height", 20)
      .append("line").attr("x1", 0).attr("x2", 50).attr("y1", 10)
      .attr("y2", 10)
      .attr("class", function (d) { return d.cls; })
      .style("stroke-width", 5);
    cells.append("p")
      .style("display", "inline-block")
      .html(function (d) { return d.name; });
      */
    return 0;
  }

  // Define fxn that runs when hit "Refresh"
  function getPars ()
  {
  var boxes = d3.selectAll("#" + args.c + " .pars");
  var vals = boxes[0].map(
      function (el) {return +d3.select(el).property("value");});
  targetMuX = vals[0]; targetMuY = vals[1]; targetSigmaX = Math.abs(vals[2]);
  targetSigmaY = Math.abs(vals[3]); targetRho = vals[4];
  //getRanges();
  generate();
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

/*****************OLD DRAWING FXN **********************/
/*
  // drawing function
  function drawDen (obj)
  {
    // set the axis domains right
    obj.y.domain(getBounds('yy'));
    obj.x.domain(getBounds('xx'));

    // draw the density lines
    var lines = obj.plot.selectAll("path")
      .data(obj.data.den);

    var line = d3.svg.line()
      .x(function(d, i) { return obj.x(d.x); })
      .y(function(d, i) { return obj.y(d.y); });

    lines.enter().append("path")
      .attr("class", function (d) { return d.name; })
      .attr("d", function (d) { return line(d.data.points); });

    // now mark the means of the densities
    var pp = posteriorMeanSD();
    var means = [mu0, xbar, pp.mean];
    var names = ['prior', 'likelihood', 'posterior'];

    var lin = obj.plot.selectAll("line")
      .data(means);
      
    var minval = obj.y.domain()[0];
    lin.enter().append("line")
      .attr("class", function (d, i) {return names[i]; })
      .attr("x1", function (d, i) { return obj.x(d); })
      .attr("x2", function (d, i) { return obj.x(d); })
      .attr("y1", function (d, i) 
          { return obj.y(minval); })
      .attr("y2", function (d, i)
          { return obj.y(obj.data.den[i].data.maxdensity); });

    var tr = obj.svg.transition().duration(500);

    tr.select(".y.axis").call(obj.yAxis);
    tr.select(".x.axis").call(obj.xAxis);

    var tr1 = obj.plot.transition().duration(500);
    tr1.selectAll("path").attr('d', function (d)
        { return line(d.data.points) } );
    tr1.selectAll("line")
      .attr("x1", function (d, i) { return obj.x(d); })
      .attr("x2", function (d, i) { return obj.x(d); })
      .attr("y2", function (d, i)
          { return obj.y(obj.data.den[i].data.maxdensity); });

    return 0;
  }
  */
  /***** NEW DRAWING FXN *********************************/
  
  
  // augmented space draw function
  function drawAug (obj)
  {
    // set the axis right
    obj.y.domain(range(obj.data.aug.yy));
    obj.x.domain(range(obj.data.aug.xx));

    // draw the target density
    var line = d3.svg.line()
      .x(function(d, i) { return obj.x(d.x); })
      .y(function(d, i) { return obj.y(d.y); });

	// obj.data.aug.points is array where element is an array corresponding to each contour line
    var lines = obj.plot.selectAll("path.contour")
      .data(obj.data.aug.points);

    lines.enter().append("path")
      .attr("class", "contour")
      .attr("d", line);

    // draw the point
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
            this._current.name = 'aug'; this._current.obj = obj; return 0; });

    // draw the path leading the point to its current position
    var lead = obj.plot.selectAll("path.lead")
      .data([obj.data.metropolis.steps]);

    lead.enter().append("path")
      .attr("class", "lead")
      .attr("d", line);

    // handle the transitions
    var tr = obj.svg.transition().duration(500);
    tr.select(".y.axis").call(obj.yAxis);
    tr.select(".x.axis").call(obj.xAxis);
    tr.select("circle").attrTween("cx", pointTweenX1)
      .attrTween("cy", pointTweenY1);
    tr.select("path.lead").delay(500).attr('d', line);
    
    tr.selectAll("path.contour").attr('d', function (d)
       { return line(d) } );
   
    lead.exit().remove();

    return 0;
  }
  
  
    function moveAlongPoints (ref, axis, field)
  {
    // determine the length of the intermediate points (the number of chunks 
    // we need to partition the time interval by is that minus one)
    var len = ref.length, tlen = len - 1, dt = 1 / tlen;
    return function (t)
    {
      var ind = d3.round(1 + (t - t % dt) / dt);
      if (ind == ref.length) ind--;
      var bound1 = dt * (ind - 1), bound2 = dt * ind;
      var beginning = ref[ind - 1], end = ref[ind];
      var changedTime = (t - bound1) / (bound2 - bound1);
      var inter = d3.interpolate(axis(beginning[field]), axis(end[field]));
      return inter(changedTime);
    };
    return 0;
  }

  function pointTweenY1 (d, i, a)
  {
    var axis = this._current.obj.y;
    return augTweenY1 (this, axis, 'y', d, i, a);
  }

  function pointTweenX1 (d, i, a)
  {
    var axis = this._current.obj.x;
    return augTweenY1 (this, axis, 'x', d, i, a);
  }

  // custom transition functions for the augmented space
  function augTweenY1 (el, axis, field, d, i, a)
  {
    // get the steps object, and the current step
    var start = el._current.d;
    el._current.d = d;
    // get reference y values (along which we will make the transition
    var refAll = data.metropolis.steps; ref = [];
    // put in the old value of the position as the starting points unless
    // its already in there, or its the first iteration/step
    var curStep = data.metropolis.currentStep;
    if (curStep == 0 || curStep == nsteps) {
      ref.push(start);
    }
    // it depends whether we are going left to right, or right to left
    // if curStep is not 7, we only need to push the last two steps
    var len = refAll.length;
    for (var ii = (curStep == 7 || curStep == -1) ? 0 : d3.max([0, len - 2]);
        ii < len; ii++) {
      var oo = refAll[ii];
      ref.push(oo);
    }
    // if we have rejected the proposal in the end, return to the initial 
    // position
    var metropolis = data.metropolis;
    if (!data.metropolis.accepted) ref.push({x : metropolis.x.val, y : metropolis.y.val});
    return moveAlongPoints(ref, axis, field);
  }
  
  // a custom transition function for the points to follow curves
  function pointTweenY (d, i, a)
  {
    // first determine what plot we're working with
    var nam = this._current.name;
    // get the steps object, and the current step
    var obj = this._current.obj, axis = obj.y;
    var start = this._current.d;
    this._current.d = d;
    // get reference y values (along which we will make the transition
    var refAll;
    if (nam == 'x') refAll = data.x.points;
    if (nam == 'y') refAll = data.y.points;
    var ref = [start];
    // it depends whether we are going left to right, or right to left
    if (d.x > start.x) {
      for (var ii = 0; ii < refAll.length; ii++) {
        var oo = refAll[ii];
        if (oo.x < d.x && oo.x > start.x) ref.push(oo);
      }
    } else {
      for (var ii = refAll.length; ii > 0; ii--) {
        var oo = refAll[ii - 1];
        if (oo.x > d.x && oo.x < start.x) ref.push(oo);
      }
    }
    ref.push(d);
    return moveAlongPoints(ref, axis, 'y');
  }
  
  
  // end of new drawing fxn

  legend();
  redraw();

  getPars();
  // register an event listener for pressing the Enter key when editing numbers
  d3.selectAll("#" + args.c + " .pars").on("keyup", function ()
      {
      if (d3.event.keyCode == 13) getPars();
    });



  return 0;
}
visual()
//addScriptToHead(args.jspath+"conrec.js",visual);

return 0;
})();
