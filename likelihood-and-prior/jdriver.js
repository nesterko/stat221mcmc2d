// a sample d3 visualization for Stat 221
// Strength of prior
// sergiy 20120208


(function() {
  var args = getURIargs();

  // put in the custom style css for the visualization
  $('head').append('<link rel="stylesheet" href="' + 
    args.jspath + '/style.css">');
        
  function visual ()
{
  // global variables for likelihood and prior means and variances
  var xbar = 0, ss = 2, mu0 = 1, ss0 = 1;
  // a global variable to tell the number of observations
  var nn = 10;

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
    drawDen(den);
    return 0;
  }

  // generate data for the plots
  function generate ()
  {
    getRanges();
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
    return 0;
  }

  // variables for resolutions and ranges of the likelihood, prior, and 
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

  // nuts and bolts functions - simulations, HMC steps, contour plots etc
  // numerics related constants
  var tobs = [-5, 5];
  var posRange = [-9, 9], posResolution = 250;
  var cauchyscale = 0.5;
  var momMean = 0, momSD = 1, momResolution = 250,
      momRange = [momMean - 3 * momSD, momMean + 3 * momSD];
  // Contour plot settings
  var nLevels = 12, cThreshold = 0.05;

  function makeDensityData (ran, resolution, densityFn, normalize)
  {
    var tmp = Object();
    tmp.xx = jStat.seq(ran[0], ran[1], resolution);
    tmp.yy = tmp.xx.map(function (e)
        { return densityFn(e); });
    var maxdensity = d3.max(tmp.yy);
    var mindensity = d3.min(tmp.yy);
    if (normalize) {
      tmp.yy = tmp.yy.map(function (e) { return e / maxdensity; });
    }
    tmp.maxdensity = maxdensity;
    tmp.mindensity = mindensity;
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
  function normal (x, mu, ss)
  {
    var fun = jStat.normal.pdf;
    return fun(x, mu, ss);
  }
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
  function makePoints (xx, yy)
  {
    return jStat.seq(0, xx.length - 1, xx.length).map (function (ii)
        {
          return { x : xx[ii], y : yy[ii] };
        });
  }
  function range (ar) { return [d3.min(ar), d3.max(ar)]; }

  // generate/initialize the data object
  var data = Object();
  generate();

  // initialize our plots and data object
  var den = init(args.c, wv, h, "Functions", 'parameter', 
      "function value", data);

  // a function to add the legend
  function legend ()
  {
    // select where we'll be adding the legend to
    var leg = d3.select("#" + args.c).append("div")
     .attr('class', "ivcont userviscontrols");
    // set the width of the conainer
    // add all fields, text inputs and the button
    var data = [
    {cls : 'prior', name : 'Prior'},
    {cls : 'likelihood', name : 'Likelihood'},
    {cls : 'posterior', name : 'Posterior'}
    ],
      databut = [
      {name : "Refresh", fn : getPars}
    ];
    var inputdata = [
    {lab : '\\\\( \\bar{x} \\\\)', cls : 'xbar', val : xbar},
    {lab : '\\\\( \\sigma \\\\)', cls : 'sigma', val : ss},
    {lab : '\\\\( n \\\\)', cls : 'nn', val : nn},
    {lab : '\\\\( \\mu_0 \\\\)', cls : 'mu0', val : mu0},
    {lab : '\\\\( \\sigma_0\\\\)', cls : 'sigma0', val : ss0}];
    var inp = leg.selectAll(".labsandpars").data(inputdata)
      .enter().append('div').attr('class', 'labsandpars');
    inp.insert('p').html(function (d) {return d.lab; });
    inp.append('input').attr('class', function (d) { return 'pars ' + d.cls})
      .attr('type', 'number')
      .attr('format', function (d, i) 
          { return i == 1? '\\d+' : '\\d'; })
      .attr('value', function (d) {return d.val;});
    //var buttonsandmore = leg.append('div').attr('class', 'buttons');
    var buttons = leg
      .selectAll(".button").data(databut);
    buttons.enter().append("input")
      .attr("type", "submit").attr("class", "button")
      .attr("value", function (d) { return d.name; })
      .on("click", function (d) { return d.fn(); });
    var table = leg.append('div').attr('class', 'tabcont')
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
    return 0;
  }

  function getPars ()
  {
  var boxes = d3.selectAll("#" + args.c + " .pars");
  var vals = boxes[0].map(
      function (el) {return +d3.select(el).property("value");});
  xbar = vals[0]; ss = Math.abs(vals[1]); nn = Math.abs(vals[2]);
  mu0 = vals[3]; ss0 = Math.abs(vals[4]);
  getRanges();
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
    //tr1.selectAll("path").attr('d', function (d)
     //   { return line(d.data.points) } );
    tr1.selectAll("line")
      .attr("x1", function (d, i) { return obj.x(d); })
      .attr("x2", function (d, i) { return obj.x(d); })
      .attr("y2", function (d, i)
          { return obj.y(obj.data.den[i].data.maxdensity); });

    return 0;
  }

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
  visual();

return 0;
})();
