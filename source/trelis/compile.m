function compile()

cd 'src';

%% compile assembly interface
mex inp_read_asm.cpp inp_assembly_section.cpp

%% compile coordinates interface
mex inp_read_coords.cpp inp_nodes_section.cpp

%% compile elements interface
mex inp_read_elems.cpp inp_elements_section.cpp

%% compile nodesets interface
mex inp_read_nsets.cpp inp_nodesets_section.cpp

%% compile  properties interface
mex inp_read_props.cpp inp_properties_section.cpp

%% compile sidesets interface
mex inp_read_ssets.cpp inp_sidesets_section.cpp

%% compile read whole file interface
mex inp_read_file.cpp inp_nodes_section.cpp inp_elements_section.cpp inp_nodesets_section.cpp inp_sidesets_section.cpp inp_properties_section.cpp inp_assembly_section.cpp

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
