// This is the builder.sce
// must be run from this directory

// [1] generate Path.incl

ilib_path_incl()

// [2] the source 

if c_link('libdistmesh_Interf') then
  printf('please do not rebuild a shared library while it is linked\n')
  printf('use ulink to unlink first\n');
else 
  printf('building shared library\n')
  chdir('src');
  // we need to chdir to execute a builder.
  ok=exec('builder.sce',errcatch=%t);
  if ~ok then 
    x_message('Compilation of source file failed\n");
  end
  chdir("../");
end 

// macros 
//--------

add_lib('macros',compile=%t);
