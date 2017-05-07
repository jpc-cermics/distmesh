libname  = 'libnspqhull' 	
libtitle='distmesh:';

// macros. 

add_lib('macros',compile=%t);

// loader for src 
exec('src/loader.sce');

printf(libtitle+' loaded\n');

// path to here 
// qhullnsp_path = get_current_exec_dir()
