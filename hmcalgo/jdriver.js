// a sample d3 visualization for Stat 221
// sergiy 20121125


(function() {
  var args = getURIargs();

  // put in the custom style css for the visualization
  $('head').append('<link rel="stylesheet" href="' + 
    args.jspath + '/style.css">');
        
  function visual ()
{

  // compute the total dimensions of the plot as provided by the 
  // URI arguments (computeDims lives in the *utils.js file)
  var dims = computeDims(args.c, args.w, args.h);

  var w = dims.w,
h = dims.h;

// define plot margins for the plots (same for all plots)
  var m = {t : 10, r : 10, b : 40, l : 50};

  // define individual svg container width, and also the width and height
  // of the plotting g elements within each svg
  var wv = w * 0.5 - 10, pw = wv - m.r - m.l, ph = h - m.t - m.b;

  function init (container, w, h, title, xlab, ylab, data)
  {
    // initialize the SVG container
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
        xAxis = d3.svg.axis().scale(x).tickSize(3, 3, 1)
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
    drawPos(pos);
    drawMom(mom);
    drawAug(aug);
    return 0;
  }

  // gnerate data for the plots
  function generate ()
  {
    // calculate target density data
    data.pos = makeDensityData(posRange, posResolution, targetDensity, true);
    data.mom = makeDensityData(momRange, momResolution, momentumDensity, false);
    data.aug = makeContourData(data.pos, data.mom, nLevels, cThreshold);

    // initialize HMC
    hmcInit();

    // calculate momentum density data
    return 0;
  }

  // nuts and bolts functions - simulations, HMC steps, contour plots etc
  // numerics related constants
  var tobs = [-4, 4];
  var posRange = [-8, 8], posResolution = 250;
  var cauchyscale = 1;
  var momMean = 0, momSD = 1, momResolution = 250,
      momRange = [momMean - 3 * momSD, momMean + 3 * momSD];
  // Contour plot settings
  var nLevels = 12, cThreshold = 0.05;

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

  function makeDensityData (ran, resolution, densityFn, normalize)
  {
    var res = Object();
    res.xx = jStat.seq(ran[0], ran[1], resolution);
    res.yy = res.xx.map(function (e)
        { return densityFn(e); });
    var maxdensity = d3.max(res.yy);
    var mindensity = d3.min(res.yy);
    if (normalize) {
      res.yy = res.yy.map(function (e) { return e / maxdensity; });
    }
    res.maxdensity = maxdensity;
    res.mindensity = mindensity;
    res.points = makePoints(res.xx, res.yy);
    return res;
  }
  function targetDensity (par)
  {
    // calculate the (unnormalized) target density as a likelihood function of 
    // a Cauchy based on 2 observations
    var pdf = jStat.cauchy.pdf;
    var vals = tobs.map(function (e) { return pdf(e, par, cauchyscale); });
    return jStat.product(vals);
  }
  function logTargetDensity (par) { return Math.log(targetDensity(par)); }
  function logTargetGradient (par)
  {
    var vals = tobs.map(function (e)
        {
          return (e - par) / (Math.pow(e - par, 2) + Math.pow(cauchyscale, 2));
        });
    return d3.sum(vals);
  }
  function momentumDensity (mo)
  {
    var pdf = jStat.normal.pdf;
    return pdf(mo, momMean, momSD);
  }
  function makePoints (xx, yy)
  {
    return jStat.seq(0, xx.length - 1, xx.length).map (function (ii)
        {
          return { x : xx[ii], y : yy[ii] };
        });
  }
  function range (ar) { return [d3.min(ar), d3.max(ar)]; }

  // HMC portion of the code
  var nsteps = 7, stepsize = 0.3;
  function hmcInit ()
  {
    data.hmc = Object();
    data.hmc.pos = Object();
    data.hmc.mom = Object();
    data.hmc.pos.val = data.hmc.mom.val = 0;
    data.hmc.steps = [ {x : 0, y : 0} ];
    data.hmc.pos.vals = [];
    data.hmc.mom.vals = [];
    data.hmc.currentStep = -1;
    data.hmc.currentIteration = 0;
    return 0;
  }
  function hmcU (par) { return -1 * logTargetDensity(par); }
  function hmcGradU (par) { return -1 * logTargetGradient(par); }
  function hmcNextStep () { var r = hmcJump (false); redraw(); }
  function hmcNextIteration () { var r = hmcJump (true); redraw(); }
  function hmc100Iterations ()
  {
    for (var ii = 0; ii < 100; ii++ ) hmcJump (true);
    redraw();
  }
  function hmcReset () { hmcInit(); redraw();}
  function hmcJump (nextIteration)
  {
    var curStep = data.hmc.currentStep;
    // declare a variable for momentum
    var pp = data.hmc.mom.val;
    // current position
    var qq = data.hmc.pos.val;
    // generate new momentum in case we are on step -1, or done with the 
    // previous iteration
    if (curStep == -1 || curStep == nsteps) {
      curStep = -1;
      pp = jStat.normal.sample(momMean, momSD);
      curStep++;
      data.hmc.currentStep = curStep;
      data.hmc.mom.val = pp;
      // initialize the object for the small history of steps
      data.hmc.steps = [{x : qq, y : pp}];
      if (!nextIteration) return 0;
    }
    // determine the number of steps we need to simulate now (go all the way 
    // to the end of the sequence if we need to go to next iteration)
    var num = 1 + nextIteration * (nsteps - curStep - 1);
    for (var ii = 0; ii < num; ii++) {
      // on the first step, we make a half-step for the momentum
      if (curStep == 0 && ii == 0) {
        pp = (pp - 0.5 * stepsize * hmcGradU(qq));
      }
      // update position with a fulls step
      qq = qq + 1 / Math.pow(momSD, 2) * pp * stepsize;

      // update momentum unless this is the last iteration, and a half step
      // in this case
      if (curStep != nsteps) {
        pp = pp - stepsize * hmcGradU(qq);
      } else {
        pp = pp - 0.5 * stepsize * hmcGradU(qq);
        // negate the momemtum at the last iteration to make the proposal 
        // symmetric
        // pp = -1 * p
      }

      // do the bookkeeping
      data.hmc.pos.val = qq;
      data.hmc.mom.val = pp;
      curStep++;
      data.hmc.currentStep = curStep;
      data.hmc.steps.push({x : qq, y : pp});
    }
    // if the current step is equal to nsteps, perform acceptance/rejection 
    if (curStep == nsteps) {
      var currentU = hmcU(data.hmc.steps[0].x),
        currentK = (0.5 * Math.pow(data.hmc.steps[0].y, 2)
            / Math.pow(momSD, 2));
      var proposedU = hmcU(qq),
          proposedK = (0.5 * Math.pow(pp, 2) / Math.pow(momSD, 2));
      var logunif = Math.log(jStat.uniform.sample(0, 1));
      if (logunif < currentU - proposedU + currentK - proposedK) {
        // accept
        data.hmc.accepted = true;
        data.hmc.pos.val = qq;
        data.hmc.pos.vals.push(qq);
        data.hmc.mom.vals.push(pp);
      } else {
        // reject
        data.hmc.accepted = false;
        data.hmc.pos.val = data.hmc.steps[0].x;
        data.hmc.pos.vals.push(data.hmc.steps[0].x);
        data.hmc.mom.vals.push(data.hmc.steps[0].y);
      }
      data.hmc.currentIteration++;
    }
    return 0;
  }
  // end of nuts and bolts functions

  // initialize our plots and data object
  var data = Object();
  var pos = init(args.c, wv, h, "Marginal position", 'position', 
      "unnormalized density", data),
      mom = init(args.c, wv, h, "Marginal momentum", 'momentum', "density",
          data),
      aug = init(args.c, wv, h, "Joint Position and Momentum",
          'position', 'momentum', data);

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

    var lines = obj.plot.selectAll("path.contour")
      .data(obj.data.aug.points);

    lines.enter().append("path")
      .attr("class", "contour")
      .attr("d", line);

    // draw the point
    var poi = {x : obj.data.hmc.pos.val, y : obj.data.hmc.mom.val};

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
      .data([obj.data.hmc.steps]);

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
    var refAll = data.hmc.steps; ref = [];
    // put in the old value of the position as the starting points unless
    // its already in there, or its the first iteration/step
    var curStep = data.hmc.currentStep;
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
    var hmc = data.hmc;
    if (!data.hmc.accepted) ref.push({x : hmc.pos.val, y : hmc.mom.val});
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
    if (nam == 'mom') refAll = data.mom.points;
    if (nam == 'pos') refAll = data.pos.points;
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

  // momentum view draw function
  function drawMom (obj)
  {
    return drawHelperPosMom('mom', obj);
  }

  function drawHelperPosMom (nam, obj)
  {
    // set the axis right
    obj.y.domain(range(obj.data[nam].yy));
    obj.x.domain(range(obj.data[nam].xx));

    // determine the point data
    var poi, redline = [0, 0];
    if (nam == 'mom') {
      poi = {x : obj.data.hmc.mom.val,
        y : momentumDensity(obj.data.hmc.mom.val)};
      redline[0] = obj.data.mom.mindensity;
      redline[1] = obj.data.mom.maxdensity;
    } else if (nam == 'pos') {
      poi = {x : obj.data.hmc.pos.val,
        y : targetDensity(obj.data.hmc.pos.val) / obj.data.pos.maxdensity};
      redline[0] = obj.data.pos.mindensity / obj.data.pos.maxdensity;
      redline[1] = 1;
    }

    // determine the histogram data, and
    // draw the histograms of the values drawn so far. Draw them first so that the rest 
    // of the illustration is above the histograms
    var ran = (nam == 'pos') ? posRange : momRange;
    var hist = d3.layout.histogram().range(ran).bins(10);
    var histdata = hist(data.hmc[nam].vals);
    // go through the histdata object and renormalize the counts so that we 
    // can superimpose the histogram and the plots
    var maxcount = d3.max(histdata.map(function (e) { return e.y; }));
    maxcount = d3.max([maxcount, 0.01]);
    for (var ii = 0; ii < histdata.length; ii++) {
      histdata[ii].y = histdata[ii].y / maxcount * (redline[1] - redline[0]) + redline[0];
    }

    // having the histogram data, draw the histogram itself
    var rects = obj.plot.selectAll("rect.histogram").data(histdata);
    rects.enter().append("rect")
      .attr('class', 'histogram')
      .attr("x", function (d) { return obj.x(d.x); })
      .attr("width", function(d) { return obj.x(d.x + d.dx) - obj.x(d.x); })
      .attr("y", function(d) { return obj.y(redline[0]); })
      .attr("height", 0);

    // draw the target density
    var line = d3.svg.line()
      .x(function(d, i) { return obj.x(d.x); })
      .y(function(d, i) { return obj.y(d.y); });

    var lines = obj.plot.selectAll("path")
      .data([obj.data[nam].points]);

    lines.enter().append("path")
      .attr("class", "pdf")
      .attr("d", line);

    // draw the point
    var point = obj.plot.selectAll("circle")
      .data([poi]);

    point.enter().append("circle")
      .attr("class", "currentValue")
      .attr("cx", function (d, i) { return obj.x(d.x); })
      .attr("cy", function (d, i) { return obj.y(d.y); })
      .attr("r", 7)
      .each( function (d)
          { this._current = Object(); this._current.d = d;
            this._current.name = nam; this._current.obj = obj; return 0; });

    // draw a vertical line
    var lin = obj.plot.selectAll("line")
      .data([obj.data.hmc[nam].val]);
      
    lin.enter().append("line")
      .attr("class", "currentValue")
      .attr("x1", function (d, i) { return obj.x(d); })
      .attr("x2", function (d, i) { return obj.x(d); })
      .attr("y1", function (d) { return obj.y(redline[0]); })
      .attr("y2", function (d) { return obj.y(redline[1]); });

    // handle the transitions
    var tr = obj.svg.transition().duration(500);
    tr.select(".y.axis").call(obj.yAxis);
    tr.select(".x.axis").call(obj.xAxis);
    tr.select("circle").attr("cx", function (d, i) { return obj.x(d.x); })
      .attrTween("cy", pointTweenY);
    tr.selectAll('rect.histogram').attr("y", function(d) { return obj.y(d.y); })
      .attr("height", function(d) { return obj.y(redline[0]) - obj.y(d.y); });

    var tr1 = obj.plot.transition().duration(500);
    tr1.select("line").attr("x1", function (d, i) { return obj.x(d); })
      .attr("x2", function (d, i) { return obj.x(d); });

    return 0;
  }

  // position view draw function
  function drawPos (obj)
  {
    return drawHelperPosMom('pos', obj);
  }

  // a function to add the legend
  function legend ()
  {
    // select where we'll be adding the legend to
    var leg = d3.select("#" + args.c).append('div')
      .attr('class', "ivcont userviscontrols");
    // set the width of the conainer right
    leg.attr("style", "width : " + wv + "px;");
    // add all fields, text inputs and the button
    var data = [
    {cls : 'lead', name : 'HMC simulated trajectory line'},
    {cls : 'pdf', name : 'Target (position) or augmented (momentum) density function'}
    ],
      databut = [
      {name : "Next step", fn : hmcNextStep},
      {name : "Next iteration", fn : hmcNextIteration},
      {name : "100 iterations", fn : hmc100Iterations},
      {name : "Reset", fn : hmcReset}
    ];
    var table = leg.append("div").attr("class", "table");
    var trow = leg.selectAll(".trow")
      .data(data)
      .enter().append("div")
      .attr("class", "trow");
    var svgs = trow.append('div').attr('class', 'tcell').append("svg")
      .attr("width", 50).attr("height", 20);
    svgs.append("line").attr("x1", 0).attr("x2", 50).attr("y1", 10)
      .attr("y2", 10)
      .attr("class", function (d) { return d.cls; });
    var captions = trow.append("div").attr("class", "ltcell");
    captions.append("p").html(function (d) { return d.name; });
    var buttons = leg.append('div').attr('class', 'buttons')
      .selectAll(".button").data(databut);
    buttons.enter().append("input")
      .attr("type", "submit").attr("class", "button")
      .attr("value", function (d) { return d.name; })
      .on("click", function (d) { return d.fn(); });
    return 0;
  }

  generate();
  legend();
  redraw();

  return 0;
}
  addScriptToHead(args.jspath + "/conrec.js", visual);

return 0;
})();
