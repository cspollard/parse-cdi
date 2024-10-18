import numpy
from cpplot.cpplot import compare, stderr
from matplotlib import figure
from CDI import sffile

def binning(xs):
  return list(map(lambda x : x[0][0], xs)) + [ xs[-1][0][1] ]


# this is bad.
# it assumes all keys are there properly in all bins...
def tohists(xs):
  nmax = len(xs)
  allkeys = xs[0][2].keys()

  centralvalue = [ x[1][0] for x in xs ]
  allvars = { k : [ xs[i][2][k] for i in range(nmax) ] for k in allkeys }


  return centralvalue, allvars


def dictmap(f, d):
  return { k : f(x) for k , x in d.items() }


if __name__ == "__main__":
  from sys import argv

  wp = int(argv[1])
  fname =  \
    "btag_ttbarPDF_mc16ade_v1.0_21-2-93_DL1r_AntiKt4EMPFlowJets_BTagging201903_%s.txt" \
    % argv[1]

  with open(fname) as reader:
    s = reader.read()

  cdi = sffile(s, 0).value
  nominal, vars = tohists(cdi)
  bins = numpy.array(binning(cdi))

  def app(h):
    return (1 + numpy.array(h))*nominal

  nominal = numpy.array(nominal)

  dvars = dictmap(app, vars)

  plotvarnames = \
    [ "Flavor_Composition"
    , "Flavor_Response"
    , "EtaIntercalibration_Modelling"
    , "FT_EFF_ttbar_PowHW7"
    , "ttbar_mc_rad"
    , "Pileup_OffsetMu"
    , "Pileup_RhoTopology"
    , "JER_DataVsMC_MC16"
    , "JER_EffectiveNP_1"
    ]

  remove = [ "PDF" ]

  vars = [ dvars[k] for k in dvars if all(r not in k for r in remove) ]
  statvars = [ dvars[k] for k in dvars if "stat" in k ]
  plotvars = [ dvars[k] for k in plotvarnames ]
  systvars = \
    [ dvars[k] for k in dvars
        if "stat" not in k
          and "singletop" not in k
          and "JER" not in k
          and "PDF" not in k
          and "Light" not in k
    ]

  alluncerts = stderr(nominal, vars)
  withstatuncerts = stderr(nominal, statvars)
  withsysuncerts = stderr(nominal, systvars)
  otcalibuncerts = stderr(nominal, plotvars)

  from maltesfs import csvsfs
  
  # for whatever reason, numpy removes "-" from recarray key names...
  uncerts = \
    [ "Flavor_Compositiondown"
    , "Pileup_OffsetMudown"
    , "JER_EffectiveNP_1down"
    , "Flavor_Responsedown"
    , "JER_DataVsMC_MC16down"
    , "EtaIntercalibrationdown"
    , "Pileup_RhoTopologydown"
    , "20230622_145054082514_hdamp"
    , "20230828_144908518131_lights_on_nominal"
    , "20230831_113912209609_herwig_pythia_MC"
    , "stats"
    ]


  otsfs = csvsfs[wp]
  xs = otsfs["bins"]
  xerrs = numpy.stack([otsfs["width"]]*2)
  ys = otsfs["center"]

  yuncerts = { k : numpy.sqrt(otsfs[k].T) for k in uncerts }
  yerrs = numpy.sqrt(numpy.sum([ otsfs[k].T for k in uncerts ], axis=0))
  yerrs = numpy.stack([yerrs]*2)


  bins = numpy.array(bins)
  bincenters = (bins[1:] + bins[:-1]) / 2.0
  binerrs = (bins[1:] - bins[:-1]) / 2.0


  fig = figure.Figure((4.5, 4.5))
  plt = fig.add_subplot(111)

  xzeros = numpy.zeros_like(xerrs)
  yzeros = numpy.zeros_like(yerrs)
  curves = \
    [ ( ys , yerrs )
    , ( ys , numpy.stack([yuncerts["stats"]]*2) )
    # , ( ys + yuncerts["20230831_113912209609_herwig_pythia_MC"] , yzeros )
    # , ( ys + yuncerts["20230828_144908518131_lights_on_nominal"] , yzeros )
    # , ( ys + yuncerts["Flavor_Compositiondown"] , yzeros )
    # , ( ys + yuncerts["Pileup_RhoTopologydown"] , yzeros )
    # , ( ys + yuncerts["JER_EffectiveNP_1down"] , yzeros )
    ]

  curvelabels = \
    [ "optimal transport"
    , "stat uncertainty"
    # , "parton shower uncertainty"
    # , "light calib uncertainty"
    # , "JES flavor composition"
    # , "JES rho topology"
    # , "JER effective NP 1"
    ]

  stdcalibs = [ alluncerts , withstatuncerts , ]
  stdlabels = [ "standard calib" , "standard calib stat" , "standard calib JER/stop" ]
  nstd = len(stdcalibs)

  nextracurves = len(curves)-2
 
  compare \
    ( plt
    , [ ( xs , xzeros ) ] * len(curves) + [ ( bincenters , binerrs ) ] * nstd
    , curves + stdcalibs
    , curvelabels + stdlabels
    , "jet $p_\\mathrm{T}$ [GeV]"
    , "efficiency scale factor"
    , alphas = [0.5] + [ 1.0 ] * 2 + [ 0.75 ] * (nstd - 1)
    , errorfills = [ True ] * len(curves) + [ False ] * nstd
    , markers = [ None ] * len(curves) + [ "" ] * nstd
    , linewidths = [ 2 , 0 ] + [ 2 ] * nstd
    , colors= \
        [ "red" ] * 2
      + [ "blue" , "green" , "magenta" , "orange" ][:nextracurves]
      + [ "black" , "gray" ][:nstd]
    )

  plt.set_xscale("log")
  plt.set_ylim(0.75, 1.05)
  plt.text(60, 0.90, "ATLAS", weight="bold", style="italic", size="large")
  plt.text(60, 0.885, "Internal", style="italic", size="large")
  plt.text(60, 0.87, "$\\sqrt{s} = 13.6$ TeV", size="large")
  plt.legend(loc="lower center")

  fig.tight_layout()
  fig.savefig("%2dWP.pdf" % wp)
