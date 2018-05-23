function compile()

cd 'src';

mex -v GCC='/usr/bin/gcc-4.9' msh_read_file.cpp msh_nodes_section.cpp msh_elements_section.cpp

if ispc
    dest_folder='win';
elseif ismac
    dest_folder='mac';
elseif isunix
    dest_folder='unix';
end

if ~exist(dest_folder,'file')
    mkdir(dest_folder);
end

movefile('*.mex*64',fullfile('..',dest_folder));

cd '..';

end
