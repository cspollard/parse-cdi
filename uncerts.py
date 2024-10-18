import numpy
from cpplot.cpplot import compare, stderr, zeroerr
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
  xs = (bins[1:] + bins[:-1]) / 2.0
  xerrs = (bins[1:] - bins[:-1]) / 2.0

  def app(h):
    return (1 + numpy.array(h))*nominal

  nominal = numpy.array(nominal)

  dvars = dictmap(app, vars)


  goodvars = { k : dvars[k] for k in dvars if "stat" not in k and "singletop" not in k }
  systuncert = stderr(nominal, goodvars.values())[1]

  fracuncerts = { k : stderr(nominal, [goodvars[k]])[1] / systuncert for k in goodvars }

  ordks = \
    sorted \
    ( fracuncerts.keys()
    , key=lambda k: numpy.max(fracuncerts[k])
    , reverse=True
    )[:15]

  fracuncerts = [ fracuncerts[k] for k in ordks ]
  sumuncert = numpy.sqrt(sum([frac*frac for frac in fracuncerts]))


  bins = numpy.array(bins)
  bincenters = (bins[1:] + bins[:-1]) / 2.0
  binerrs = (bins[1:] - bins[:-1]) / 2.0


  fig = figure.Figure((5, 5))
  plt = fig.add_subplot(111)

  xzeros = numpy.zeros_like(xerrs)

  compare \
    ( plt
    , [ ( bincenters , binerrs ) ] * (1 + len(fracuncerts))
    , [ zeroerr(sumuncert) ] + [ zeroerr(u) for u in fracuncerts ]
    , ["quadrature sum"] + ordks
    , "jet $p_\\mathrm{T}$ [GeV]"
    , "included systematic uncertainty fraction"
    , alphas = [ 1.0 ] * (len(fracuncerts) + 1)
    , errorfills = [ True ] * (len(fracuncerts) + 1)
    , markers = [ None ] * (len(fracuncerts) + 1)
    , linewidths = [ 2 ] * (len(fracuncerts) + 1)
    , colors= ["black", "blue" , "green" , "gray" , "magenta" , "orange" ]
    )

  plt.plot([0, 1000], [1, 1], lw=1, color="black", ls="--")
  plt.set_xscale("log")
  plt.set_xlim(20, 400)
  plt.set_ylim(0, 1.5)
  plt.text(25, 1.4, "ATLAS", weight="bold", style="italic", size="large")
  plt.text(25, 1.33, "Internal", style="italic", size="large")
  plt.text(25, 1.26, "$\\sqrt{s} = 13.6$ TeV", size="large")
  plt.legend(loc="upper right")

  fig.tight_layout()
  fig.savefig("%2dWP-uncert.pdf" % wp)
