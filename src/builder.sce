// This is the builder.sce 
// must be run from this directory 

// [1] generate Path.incl 
if exists('%nsp') then 
  ilib_path_incl()
end 

// [3] the part devoted to shared lib generation 
ilib_name  = 'libdistmesh' 		// interface library name 

// objects files (but do not give mexfiles here)
files = ['dsegment.o','dellipse.o','dellipsoid.o','trisurfupd.o'];

// other libs needed for linking (must be shared library names)
libs =[];

// table of (scilab_name,interface-name or mexfile-name, type) 
table =['dsegment','int_dsegment';
	'dellipse','int_dellipse';
	'dellipsoid','int_dellipsoid';
	'trisurfupd','int_trisurfupd'];

ldflags ="-llapack"// "`pkg-config lp_solve --libs`";
cflags = ""// "`pkg-config lp_solve --cflags`"

// do not modify below 
// ----------------------------------------------

ilib_build(ilib_name,table,files,libs,ldflags = ldflags,cflags = cflags );


