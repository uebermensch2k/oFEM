function compile()

cd 'src';

%% compile assembly interface
%mex msh_read_asm.cpp msh_assembly_section.cpp

%% compile coordinates interface
%mex msh_read_coords.cpp msh_nodes_section.cpp

%% compile elements interface
%mex msh_read_elems.cpp msh_elements_section.cpp

%% compile nodesets interface
%mex msh_read_nsets.cpp msh_nodesets_section.cpp

%% compile  properties interface
%mex msh_read_props.cpp msh_properties_section.cpp

%% compile sidesets interface
%mex msh_read_ssets.cpp msh_sidesets_section.cpp

%% compile read whole file interface
mex msh_read_file.cpp msh_nodes_section.cpp msh_elements_section.cpp msh_physical_section.cpp

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
