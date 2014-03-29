

STATSFILE_Template = "${OUTDIR>/}genstats.py";

# these qualifiers are added to qualifiers passed to _writestat
STATQUALS = ()

def _clearstats ():
  if exists(STATSFILE):
    info("removing $STATSFILE");
    x.sh("rm $STATSFILE");
  
import fcntl
  
def _statsfile ():
  """Initializes stats file (if not existing), returns open file""";
  ff = file(STATSFILE,"a");
  fcntl.flock(ff,fcntl.LOCK_EX);
  # seek to end of file, if empty, make header
  ff.seek(0,2);
  if not ff.tell():
    ff.write("""# auto-generated noise stats file\n""");
    ff.write("""settings = dict(%s)\n"""%",".join([ "%s=%s"%(key,repr(val)) for key,val in globals().iteritems() 
                if not callable(val) and key.upper() == key ]));
    ff.write("noisestats = {}\n");
  return ff;
  
def _writestat (name,value,*qualifiers):
  ff = _statsfile();
  subsets = list(STATQUALS) + list(qualifiers);
  subsets = ','.join(map(repr,subsets));
  ff.write("noisestats.setdefault('%s',{})[%s] = %s\n"%(name,subsets,repr(value)));
  fcntl.flock(ff,fcntl.LOCK_UN);
