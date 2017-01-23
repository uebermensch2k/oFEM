

classdef mesh < handle
    properties (Constant=true)
        % vaccum permittivity
        eps0 = 8.854187817620e-12; % in A^2s^4/(kg m^3)
        % vacuum permeability
        mu0  = 1.2566370614e-6;    % in kg*m*s^-2*A^-2
    end

    properties (SetAccess=protected)
        % dimension of the space described by the mesh
        dim;

        % ofem.matrixarray of size(co,3) length containing the coordinates
        % of the mesh
        co;

        % ed is a matrix size(ed,2)=2 containing the edge to nodes mapping.
        % More precisely, each row contains the ids of the two nodes that
        % represent the respective edge. The ids are sorted in ascending
        % order.
        ed;

        % 
        fa;

        % el is a matrix with dim+1=size(el,2) containing the
        % element to coordinates mapping. More precisely, each row contains
        % the indices into the co matrix representing the coordinates of
        % its vertices.
        el;

        % el2ed is a matrix with size(el2ed,2) dependent on the type of
        % element. E.g., for tetrahedra meshes size(el2ed,2)=4. It contains
        % the element to edge mapping. More precisely, per row the ids of
        % the edges the respective elements is build up from are present.
        el2ed;

        % el2edsign is a matrix the same size as el2ed. Since ed is sorted
        % in ascending order in each row this matrix contains the
        % information if the direction of a particular edge of a particular
        % element is the same (value equals 1) direction or opposite (value
        % equals -1) direction as compared to ed.
        el2edsign

        % parts is a cell array with 3=size(parts,1) containing in each
        % column the name of the element set, the index into the material
        % data mat describing the material the element set is build from
        % and a vector of indices into el describing which elements
        % the element set is created from.
        parts;

        % bd is a cell array with 2=size(bd,1) containing in each column
        % the name of the set and a cell array. For dim==2 this cell array
        % is 2x3 and for dim==3 it is 2x4. In any case the first row
        % contains the names of the subset. The second row contains vectors
        % of indices of elements, i.e., indices into el. For dim==2
        % the boundary bd{2,1} is 'The boundary created taking the first
        % edge of the elements with indices bd{2,1}{2,1}, following the 
        % second edge of the elements with indices bd{2,1}{2,2} and the
        % third  edge of the elements with indices bd{2,1}{2,3}'.
        % The procedure is analogous for dim==2, however, a tetrahedron is
        % build up from 4 faces.
        bd;

        % type is a string containing one of the following values
        %   'line' : The elements are edges
        %   'tri'  : The elements are triangles
        %   'tet'  : The elements are tetrahedrons
        %   'quad' : The elements are convex quadrilateral
        %   'hex'  : The elements are convex hexahedrons
        type;

        % nodes_per_elem contains the number of nodes a single element is
        % build up from. This information is needed as, e.g., when
        % extending a mesh for use with P2 elements this information is
        % lost.
%         nodes_per_elem;

        % fe is the finite element used for discretization
        fe;
    end

    properties
        % mat is a cell array containing the materials the meshed model is
        % build up from. Each entry in the array is a structure with a
        % field name containing the name of the material. Dependent on the
        % material there can be additional fields, all uppercase with
        % whitespaces replaced by underscores, containing the desired
        % property as, e.g., DENSITY. parts{2,:} containg the indices
        % into mat.
        mat;
    end

    properties (SetAccess=protected,GetAccess=public)
        % Nco is the number of nodes in the loaded mesh. This number is not
        % changed by operations implicitly extending the co
        % ofem.matrixarray.
        Nco;
    end

    methods (Access=private)
        %%
        function [co,el]=create_midpoints(obj)
            Nel=numel(obj.el(:,1));
            
            co = obj.co;

            switch obj.type
                case 'edge'
                    %% edge
                    % create new nodes
                    co2 = (obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,2)))/2;
                    co2 = permute(co2,[3,1,2]); % coordinate per row
                    
                    % get rid of doubles
                    [co2,~,ib]=unique(co2,'rows','stable');
                    
                    % add new indices
                    Nco2    = size(co2,1);
                    idx2    = obj.Nco+(1:Nco2);
                    
                    % extend nodes
                    co(:,:,end+1:end+Nco2) = permute(co2,[2,3,1]);
                    
                    % append new nodes to elements
                    el=[obj.el, reshape(idx2(ib),Nel,[])];
                    
                case 'tri'
                    %% triangle
                    % create new nodes
                    mpt1 = permute(obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,2)),[3,1,2]);
                    mpt2 = permute(obj.co(:,:,obj.el(:,2))+obj.co(:,:,obj.el(:,3)),[3,1,2]);
                    mpt3 = permute(obj.co(:,:,obj.el(:,3))+obj.co(:,:,obj.el(:,1)),[3,1,2]);
                    
                    co2 = [ mpt1;  mpt2;  mpt3 ]/2;
                    
                    % get rid of doubles
                    [co2,~,ib]=unique(co2,'rows','stable');
                    
                    % add new indices
                    Nco2    = size(co2,1);
                    idx2=obj.Nco+(1:Nco2);
                    
                    % extend nodes
                    co(:,:,end+1:end+Nco2) = permute(co2,[2,3,1]);
                    
                    % append new nodes to elements
                    el=[obj.el, reshape(idx2(ib),Nel,[])];
                    
                case 'quad'
                    %% quadrilateral
                    error('ofem:mesh:NotImplemented',...
                          'Quadrilateral meshes not supported so far!');
                    
                case 'tet'
                    %% tetrahedron
                    % create new nodes
                    mpt1 = permute(obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,2)),[3,1,2]);
                    mpt2 = permute(obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,3)),[3,1,2]);
                    mpt3 = permute(obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,4)),[3,1,2]);
                    mpt4 = permute(obj.co(:,:,obj.el(:,2))+obj.co(:,:,obj.el(:,3)),[3,1,2]);
                    mpt5 = permute(obj.co(:,:,obj.el(:,2))+obj.co(:,:,obj.el(:,4)),[3,1,2]);
                    mpt6 = permute(obj.co(:,:,obj.el(:,3))+obj.co(:,:,obj.el(:,4)),[3,1,2]);
                    
                    co2 = [ mpt1; mpt2; mpt3; mpt4; mpt5; mpt6 ]/2;
                    
                    % get rid of doubles
                    [co2,~,ib]=unique(co2,'rows','stable');
                    
                    % add new indices
                    Nco2    = size(co2,1);
                    idx2=obj.Nco+(1:Nco2);
                    
                    % extend nodes
                    co(:,:,end+1:end+Nco2) = permute(co2,[2,3,1]);
                    
                    % append new nodes to elements
                    el=[obj.el, reshape(idx2(ib),Nel,[])];
                    
                case 'hex'
                    %% hexahedron
                    error('ofem:mesh:NotImplemented',...
                          'Hexahedral meshes not supported so far!');
                    
                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
        end
    end

    methods
        %%
        function valid=is_valid(obj)
            valid = ~isempty(obj.type);
        end

        %%
        function reset(obj)
        %RESET reset the mesh to its initial configuration.
        
            if ~obj.is_valid()
                obj.co    = [];
                obj.el    = [];
                obj.bd    = [];
                obj.dim   = [];
                obj.parts = [];
                obj.Nco   = [];
                obj.mat   = [];
            else
                obj.co = obj.co(:,:,1:obj.Nco);
                obj.el = obj.el(:,1:obj.nodes_per_element());
            end

            obj.ed = [];
            obj.fa = [];
            

            obj.el2ed = [];
            obj.el2edsign = [];

            obj.fe = [];
        
        end

        %%
        function noperel = nodes_per_element(obj)
            if ~obj.is_valid()
                error('ofem:mesh:InvalidMesh',...
                      'nodes_per_element called on invalid mesh');
            end

            switch obj.type
                case 'edge'
                    %% edge
                    noperel = 2;

                case 'tri'
                    %% triangle
                    noperel = 3;

                case 'quad'
                    %% quadrilateral
                    noperel = 4;

                case 'tet'
                    %% tetrahedron
                    noperel = 4;

                case 'hex'
                    %% hexahedron
                    noperel = 8;

                otherwise
                    error('ofem:mesh:InvalidMesh',...
                          'nodes_per_element called on invalid mesh');
            end
        end

        %%
        function hypercube(obj,x0,varargin)
        %hypercube creates a hypercube
        %
        %  hypercube(x0) creates a unit hypercube with lower-left corner at
        %  x0. The dimension of the embedding space is numel(x0).
        %
        %  hypercube(x0,extent) additionally sets the extent of the
        %  hypercube. extent can either be a scalar, in which case the
        %  hypercube has all equal lengths, or must match the dimension of
        %  x0.
        %
        %  In 2D, the mesh consists of two triangles, in 3D it consists of
        %  6 tetrahedra.
        %
        
            p=inputParser;
            addRequired(p,'x0',@(x) validateattributes(x,{'numeric'},{'nonempty','real','vector'}));
            addOptional(p,'extent',1,@(x) validateattributes(x,{'numeric'},{'nonempty','real','vector'}));
            parse(p,x0,varargin{:});

            extent=p.Results.extent;

            if any(extent<=0)
                error('ofem:mesh:create:argchk','The extent must be positive');
            end

            if numel(extent)==1
                extent=repmat(extent,numel(x0),1);
            else
                if numel(x0)~=numel(extent)
                    error('ofem:mesh:create:argchk','numel(extent) must match numel(x0) or be scalar');
                end
            end

            LL=x0;
            UR=x0+extent;

            obj.reset();

            obj.dim = numel(LL);
            
            switch obj.dim
                case 1
                    %% coordinates
                    obj.Nco = 2;
                    obj.co = [LL(1); UR(1)];

                    %% elements
                    obj.el = [ 1,2 ];                    

                    %% element type
                    obj.type = 'edge';
                    
                    %% sidesets
                    obj.bd = {'boundary';{'boundary_P1','boundary_P2'; 1,1}};
                    
                case 2
                    %% coordinates
                    obj.Nco = 4;
                    obj.co  = [        ...
                        [LL(1), LL(2)]; ...
                        [UR(1), LL(2)]; ...
                        [UR(1), UR(2)]; ...
                        [LL(1), UR(2)]  ...
                        ];
                    
                    %% elements
                    obj.el = [ ...
                        1,2,4; ...
                        2,3,4  ...
                        ];                    

                    %% element type
                    obj.type = 'tri';
                    
                    %% sidesets
                    obj.bd = {'boundary';{'boundary_E1','boundary_E2','boundary_E3'; [1;2],2,1}};
                    
                case 3
                    %% coordinates
                    obj.Nco = 8;
                    obj.co  = [        ...
                        [LL(1), LL(2), LL(3)]; ...
                        [UR(1), LL(2), LL(3)]; ...
                        [UR(1), UR(2), LL(3)]; ...
                        [LL(1), UR(2), LL(3)]; ...
                        [LL(1), LL(2), UR(3)]; ...
                        [UR(1), LL(2), UR(3)]; ...
                        [UR(1), UR(2), UR(3)]; ...
                        [LL(1), UR(2), UR(3)]; ...
                        ];
                    
                    %% elements
                    obj.el = [ ...
                        1,2,3,7; ...
                        1,3,4,7; ...
                        1,2,7,6; ...
                        5,1,8,7; ...
                        1,4,8,7; ...
                        5,1,7,6; ...
                        ];
                    
                    %% element type
                    obj.type = 'tet';

                    %% sidesets
                    % S1: 1, 2, 3
                    % S2: 1, 2, 4
                    % S3: 2, 3, 4
                    % S4: 1, 3, 4
                    obj.bd = {'boundary';{'boundary_S1','boundary_S2','boundary_S3','boundary_S4'; [1;2;4;5],[3;6],[1;2;3;5],[4;6]}};
                    
                otherwise
                    error('ofem:mesh:InvalidMesh',...
                          'hypercube requested for dim>3');
            end

            obj.co = ofem.matrixarray(permute(obj.co,[2,3,1]));

            %% parts
            obj.parts = {'hypercube';1;(1:size(obj.el,1))'};

            %% materials
            obj.mat{1}.name = 'air';
            obj.mat{1}.ELASTIC = [0;0];
            obj.mat{1}.DENSITY = 0.001;
            obj.mat{1}.CONDUCTIVITY = 0;
            obj.mat{1}.SPECIFIC_HEAT = 1.012;
        end

        %%
        function uniform_refine(obj)
        %uniform_refine uniformly refines the mesh
        %
            if ~obj.is_valid()
                error('ofem:mesh:InvalidMesh','ofem.mesh needs to hold a valid triangulation');
            end

            recover_fe=0;

            if ~isempty(obj.fe)
                fe_tmp=obj.fe;
                obj.reset();
                recover_fe=1;
            end
            

            [obj.co,elem] = obj.create_midpoints();

            obj.Nco = size(obj.co,3);

            Nel = size(obj.el,1);

            switch obj.type
                case 'edge'
                    %% edge
                    % update elements
                    obj.el = [         ...
                        elem(:,[1,3]); ...
                        elem(:,[3,2])  ...
                        ];

                    % update parts
                    obj.parts(3,:)=cellfun(@(x) [x;x+Nel],obj.parts(3,:),'UniformOutput',false);

                    % update boundaries
                    obj.bd{2,:}(2,2)=cellfun(@(x) x+Nel,obj.bd{2,:}(2,2),'UniformOutput',false);

                case 'tri'
                    %% triangle
                    % update elements
                    obj.el = [ ...
                        elem(:,[1,4,6]); ...
                        elem(:,[4,2,5]); ...
                        elem(:,[6,5,3]); ...
                        elem(:,[4,5,6])  ...
                        ];

                    % update parts
                    obj.parts(3,:)=cellfun(@(x) [x;x+Nel;x+2*Nel;x+3*Nel],obj.parts(3,:),'UniformOutput',false);

                    % update boundaries
                    obj.bd{2,:}(2,1)=cellfun(@(x) [x      ; x+  Nel],obj.bd{2,:}(2,1),'UniformOutput',false);
                    obj.bd{2,:}(2,2)=cellfun(@(x) [x+  Nel; x+2*Nel],obj.bd{2,:}(2,2),'UniformOutput',false);
                    obj.bd{2,:}(2,3)=cellfun(@(x) [x+2*Nel; x      ],obj.bd{2,:}(2,3),'UniformOutput',false);

                case 'quad'
                    %% quadrilateral
                    error('ofem:mesh:NotImplemented',...
                          'Quadrilateral meshes not supported so far!');

                case 'tet'
                    %% tetrahedron
%                     error('ofem:mesh:NotImplemented',...
%                           'Tetrahedral meshes not supported so far!');

                    obj.el = [ ...
                        elem(:,[1,5, 6, 7]); ...
                        elem(:,[5,2, 8, 9]); ...
                        elem(:,[6,8, 3,10]); ...
                        elem(:,[7,9,10, 4]); ...
                        elem(:,[7,5, 6, 9]); ...
                        elem(:,[5,8, 6,10]); ...
                        elem(:,[6,8,10, 9]); ...
                        elem(:,[7,9, 6,10]); ...
                        ];

                    % update parts
                    obj.parts(3,:)=cellfun(@(x) [x;x+Nel;x+2*Nel;x+3*Nel;x+4*Nel;x+5*Nel;x+6*Nel;x+7*Nel],obj.parts(3,:),'UniformOutput',false);

                    % update boundaries
                    % S1: 1, 2, 3
                    % S2: 1, 2, 4
                    % S3: 2, 3, 4
                    % S4: 1, 3, 4
                    obj.bd{2,:}(2,1)=cellfun(@(x) [x      ; x+1*Nel; x+2*Nel; x+5*Nel],obj.bd{2,:}(2,1),'UniformOutput',false);
                    obj.bd{2,:}(2,2)=cellfun(@(x) [x      ; x+1*Nel; x+3*Nel; x+4*Nel],obj.bd{2,:}(2,2),'UniformOutput',false);
                    obj.bd{2,:}(2,3)=cellfun(@(x) [x+1*Nel; x+2*Nel; x+3*Nel; x+6*Nel],obj.bd{2,:}(2,3),'UniformOutput',false);
                    obj.bd{2,:}(2,4)=cellfun(@(x) [x      ; x+2*Nel; x+3*Nel; x+7*Nel],obj.bd{2,:}(2,4),'UniformOutput',false);

                case 'hex'
                    %% hexahedron
                    error('ofem:mesh:NotImplemented',...
                          'Hexahedral meshes not supported so far!');

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end

            if recover_fe==1
                obj.assign_fe(fe_tmp);
            end
        end

        %%
        function info=load_from_inp(obj,inp_file_name)
        %LOAD_FROM_INP loads a inp-file into the ofem.mesh class.
        %
        % LOAD_FROM_INP(inp_file_name) loads the inp-file specified by
        % inp_file_name into the ofem.mesh class.
        %

            % check if extension was specified
            [filepath,filename,fileext] = fileparts(inp_file_name);
            if isempty(fileext)
                inp_file_name = fullfile(filepath,strcat(filename,'.inp'));
            end
            if isempty(filepath)
                inp_file_name = fullfile(pwd,inp_file_name);
            end

            % check if inp file exists
            if ~exist(inp_file_name,'file') % for whatever reason this fails on the cluster
                error('ofem:mesh:FileNotFound','inp-file "%s" not found.\n',inp_file_name);
            end

            tic
            inp=inp_read_file(inp_file_name);
            info.time.file_read=toc;

            % inp is a cell array of the following structure
            % inp{:,1}: coordinates
            % inp{:,2}: elements
            % inp{:,3}: nodesets
            % inp{:,4}: sidesets
            % inp{:,5}: properties
            % inp{:,6}: materials

            tic

            %% coordinates
            % size(obj.co)=[Nd,1,Nco], always column vectors
            obj.co  = ofem.matrixarray(permute(inp{2,1}{2,1},[2,3,1]));
            obj.Nco = size(obj.co,3);
            obj.dim = size(obj.co,1);

            %% elements
            elaux  = inp{2,2};
            obj.el = vertcat(elaux{2,:});

            %% parts
            Nparts    = numel(elaux(1,:));
            obj.parts = cell(3,Nparts);

            % set name
            obj.parts(1,:) = elaux(1,:);
            % and indices
            max_idx=0;
            for i=1:Nparts
                Nel            = size(elaux{2,i},1);
                obj.parts{3,i} = max_idx+(1:Nel)';
                max_idx        = max_idx+Nel;
            end
            clear elaux;

            %% element type
            nodes_per_elem=size(obj.el,2);
            switch obj.dim
                case 1
                    if nodes_per_elem~=2
                        error('ofem:mesh:InvalidMesh',...
                              'In 1D the elements are expected to contain two nodes, but %d nodes were found.',...
                              nodes_per_elem);
                    end
                    obj.type = 'line';

                case 2
                    switch nodes_per_elem
                        case 3
                            obj.type = 'tri';
                        case 4
                            obj.type = 'quad';
                        otherwise
                            error('ofem:mesh:InvalidMesh',...
                              'In 2D the elements are expected to contain three of four nodes, but %d nodes were found.',...
                              nodes_per_elem);
                    end

                case 3
                    switch nodes_per_elem
                        case 4
                            obj.type = 'tet';
                        case 8
                            obj.type = 'hex';
                        otherwise
                            error('ofem:mesh:InvalidMesh',...
                              'In 3D the elements are expected to contain four of eight nodes, but %d nodes were found.',...
                              nodes_per_elem);
                    end
            end

            %% nodesets
%             bd_ns_name=inp{2,3};

            %% sidesets
            obj.bd=inp{2,4};

            %% materials
            mataux=inp{2,6};
            if isempty(mataux)
              % if no material present => set air as default material
              obj.mat{1}.name = 'air';
              obj.mat{1}.ELASTIC = [0;0];
              obj.mat{1}.DENSITY = 0.001;
              obj.mat{1}.CONDUCTIVITY = 0;
              obj.mat{1}.SPECIFIC_HEAT = 1.012;

              for i=1:Nparts
                obj.parts{2,i}=1;
              end
            else
              % else read materials
              Nmat=size(mataux,2);
              obj.mat=cell(1,Nmat);
              
              for i=1:Nmat
                obj.mat{i}.name=mataux{1,i};
                Nmatpar=numel(mataux{2,i}(1,:));
                for j=1:Nmatpar
                  obj.mat{i}.(strrep(mataux{2,i}{1,j},' ','_'))=mataux{2,i}{2,j};
                end
                
              end
              clear mataux;
              
              % and associate them with the element sets
              props=inp{2,5};
              Nprops=numel(props(1,:));
              
              matnames = cellfun(@(c) c.name,obj.mat,'UniformOutput',0);
              elnames  = obj.parts(1,:);
              
              for i=1:Nprops
                elidx  = strcmp(props{2,i}{1},elnames );
                matidx = find(strcmp(props{2,i}{2},matnames));
                
                obj.parts{2,elidx}=matidx;
              end
            end
            
            info.time.post_proccess=toc;

            clear inp;
        end
        
        %%
        function create_edges(obj)
        %CREATE_EDGES creates edges information and element to edges mapping
        %
            Nt = size(obj.el,1);

            switch obj.type
                case 'edge'
                    %% edge
                    obj.ed        = obj.el;
                    obj.el2ed     = 1:Nt;
                    obj.el2edsign = ones(Nt,1);
                    return;

                case 'tri'
                    %% triangle
                    obj.ed = [obj.el(:,[1 2]); ...
                              obj.el(:,[2 3]); ...
                              obj.el(:,[3 1])];

                case 'quad'
                    %% quadrilateral
                    obj.ed = [obj.el(:,[1 2]); ...
                              obj.el(:,[2 3]); ...
                              obj.el(:,[3 4]); ...
                              obj.el(:,[4 1])];

                case 'tet'
                    %% tetrahedron
                    obj.ed = [obj.el(:,[1 2]); ...
                              obj.el(:,[1 3]); ...
                              obj.el(:,[1 4]); ...
                              obj.el(:,[2 3]); ...
                              obj.el(:,[2 4]); ...
                              obj.el(:,[3 4])];

                case 'hex'
                    %% hexahedron
                    obj.ed = [obj.el(:,[1 2]); ...
                              obj.el(:,[2 3]); ...
                              obj.el(:,[3 4]); ...
                              obj.el(:,[4 1]); ...
                              obj.el(:,[5 6]); ...
                              obj.el(:,[6 7]); ...
                              obj.el(:,[7 8]); ...
                              obj.el(:,[8 5]); ...
                              obj.el(:,[1 5]); ...
                              obj.el(:,[2 6]); ...
                              obj.el(:,[3 7]); ...
                              obj.el(:,[4 8])];

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end

            [obj.ed,I]    = sort(obj.ed,2);
            I             = diff(I,1,2);%(:,2)-I(:,1);
            [obj.ed,~,ic] = unique(obj.ed,'rows','stable');
            obj.el2ed     = reshape(ic,Nt,[]);
            obj.el2edsign = int8(reshape(I,Nt,[]));
        end

        %%
        function create_faces(obj)
        %CREATE_FACES creates faces information and element to faces mapping
        %
            Nt = size(obj.el,1);

            switch obj.type
                case 'edge'
                    %% edge
                    warning('ofem:mesh:Invalid',...
                            'create_faces called but mesh is composed of edges!');
                    return;

                case 'tri'
                    %% triangle
                    obj.fa        = obj.el;
                    obj.el2fa     = 1:Nt;
                    obj.el2fasign = ones(Nt,1);

                case 'quad'
                    %% quadrilateral
                    obj.fa        = obj.el;
                    obj.el2fa     = 1:Nt;
%                     obj.el2fasign = ones(Nt,1);

                case 'tet'
                    %% tetrahedron
                    % S1: 1, 2, 3
                    % S2: 1, 2, 4
                    % S3: 2, 3, 4
                    % S4: 1, 3, 4
                    obj.fa = [obj.el(:,[1 2 3]); ...
                              obj.el(:,[1 2 4]); ...
                              obj.el(:,[2 3 4]); ...
                              obj.el(:,[1 3 4]);];

                case 'hex'
                    %% hexahedron
                    error('ofem:mesh:InvalidMesh',...
                          'Hexahedral meshes not supported so far!');

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end

            [obj.fa,I]    = sort(obj.fa,2);
            I             = I(:,2)-I(:,1);
            [obj.ed,~,ic] = unique(obj.ed,'rows','stable');
            obj.el2ed     = reshape(ic,Nt,[]);
            obj.el2edsign = INT8(reshape(I,Nt,[]));
        end

        %%
        function assign_fe(obj,fe)
            if ~obj.is_valid()
                error('ofem:mesh:InvalidMesh','ofem.mesh needs to hold a valid triangulation');
            end

            switch fe
                case ofem.finiteelement.P1
                    % nothing to be done
                    
                case ofem.finiteelement.P2
                    % extend mesh
                    [obj.co,obj.el] = obj.create_midpoints();

                case ofem.finiteelement.Q1
                    error('ofem:mesh:InvalidMesh',...
                          'Q1 elements not supported so far!');

                case ofem.finiteelement.Q2
                    error('ofem:mesh:InvalidMesh',...
                          'Q2 elements not supported so far!');

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end

            obj.fe=fe;
        end

        %%
        function [DinvT,detD]=jacobiandata(obj)
        %JACOBIANDATA returns the Jacobians and their determinants
        %
        % [DinvT,detD]=JACOBIANDATA, for a transformation
        % x=\Phi_k(y)=D_ky+n_k
        % from a reference element given in y-coordinates to the global
        % element \Omega_k given in x-coordinates, returns the Jacobians
        % D_k^{-T} and their determinants. DinvT is a matrixarray build up
        % as follows
        %   DinvT(:,:,i) = D_i^{-T}
        % likewise
        %   detD(:,:,i) = det(DinvT(:,:,i))
        %
            switch obj.type
                case 'edge'
                    %% edge
                    tmp=obj.co(:,1,obj.el(:,2))-obj.co(:,1,obj.el(:,1));
                    DinvT = 1/tmp;
                    detD  = tmp;

                case 'tri'
                    %% triangle
                    e12 = obj.co(:,1,obj.el(:,2))-obj.co(:,1,obj.el(:,1));
                    e13 = obj.co(:,1,obj.el(:,3))-obj.co(:,1,obj.el(:,1));

                    detD = dot(rot(e12),e13,1);

                    DinvT = [ -rot(e13), rot(e12) ]./detD;

                case 'quad'
                    %% quadrilateral
                    error('ofem:mesh:InvalidMesh',...
                          'Quadrilateral meshes not supported so far!');

                case 'tet'
                    %% tetrahedron
                    e12 = obj.co(:,1,obj.el(:,2))-obj.co(:,1,obj.el(:,1));
                    e13 = obj.co(:,1,obj.el(:,3))-obj.co(:,1,obj.el(:,1));
                    e14 = obj.co(:,1,obj.el(:,4))-obj.co(:,1,obj.el(:,1));
                    
                    detD = dot(e12,cross(e13,e14,1),1);

                    DinvT = [ cross(e13,e14,1) , ...
                              cross(e14,e12,1) , ...
                              cross(e12,e13,1) ]./detD;

                case 'hex'
                    %% hexahedron
                    error('ofem:mesh:InvalidMesh',...
                          'Hexahedral meshes not supported so far!');

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
        end

        %% get the Dirichlet boundary nodes
        function nID=dirichlet(obj,idx)
        %DIRICHLET returns node indices for the specified sidesets.
        %
        % nID=DIRICHLET(idx), where idx is an index vector indexing the bd 
        % cell, returns a cell array the same length as idx. In each cell
        % the node IDs of the respective sideset are present.
        %
            ss    = obj.bd(2,idx);
            Nss   = size(ss,2);
            nID   = cell(Nss,1);
            npere = size(obj.el,2);
            switch obj.type
                case 'edge'
                    %% edges

                case 'tri'
                    %% triangles
                    % sideset coding for 2D triangle structures by the
                    % Abaqus inp file is as follows:
                    % E1: 1,2
                    % E2: 2,3
                    % E3: 3,1
                    switch npere
                        case 6
                            %% second order mesh
                            % 1 2 3 e12 e23 e31
                            for i=1:Nss
                                nID{i}=unique([obj.el(ss{1,i}{2,1},[1,4,2]); ...
                                               obj.el(ss{1,i}{2,2},[2,5,3]); ...
                                               obj.el(ss{1,i}{2,3},[3,6,1])],'stable');
                            end
                        case 3
                            %% first order mesh
                            % 1 2 3
                            for i=1:Nss
                                nID{i}=unique([obj.el(ss{1,i}{2,1},[1,2]); ...
                                               obj.el(ss{1,i}{2,2},[2,3]);...
                                               obj.el(ss{1,i}{2,3},[3,1])],'stable');
                            end
                        otherwise
                            error('ofem:mesh:Unspecified',...
                                  'Triangle expected but the number of elemnts does not fit!\nThe number found is %d',...
                                  npere);
                    end
                            
                case 'quad'
                    %% quadrilateral
                    error('ofem:mesh:InvalidMesh',...
                        'Quadrilateral meshes not supported so far!');
                case 'tet'
                    %% tets
                    % sideset coding for 3D tetrahedral structures by
                    % the Abaqus inp file is as follows:
                    % S1: 1, 2, 3
                    % S2: 1, 2, 4
                    % S3: 2, 3, 4
                    % S4: 1, 3, 4
                    switch npere
                        case 10
                            %% second order mesh
                            % 1 2 3 4 e12 e13 e14 e23 e24 e34
                            for i=1:Nss
                                nID{i}=unique([obj.el(ss{1,i}{2,1},[1,5,2, 8,3,6]); ...
                                               obj.el(ss{1,i}{2,2},[1,5,2, 9,4,7]);...
                                               obj.el(ss{1,i}{2,3},[2,8,3,10,4,9]);
                                               obj.el(ss{1,i}{2,4},[1,6,3,10,4,7])],'stable');
                            end
                        case 4
                            %% first order mesh
                            % 1 2 3 4
                            for i=1:Nss
                                nID{i}=unique([obj.el(ss{1,i}{2,1},[1,2,3]); ...
                                               obj.el(ss{1,i}{2,2},[1,2,4]);...
                                               obj.el(ss{1,i}{2,3},[2,3,4]);
                                               obj.el(ss{1,i}{2,4},[1,3,4])],'stable');
                            end
                        otherwise
                            error('ofem:mesh:Unspecified',...
                                  'Tetrahedron expected but the number of elemnts does not fit!\nThe number found is %d',...
                                npere);
                    end
 
                case 'hex'
                    %% hexahedron
                    error('ofem:mesh:InvalidMesh',...
                          'Hexahedral meshes not supported so far!');
                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
        end

        %% get the Neumann boundary nodes
        function [meas,faces,normals,opp]=neumann(obj,idx)
        %NEUMANN returns the normals associated to the Neumann boundary.
        %
        % data=NEUMANN(idx), where idx is an index vector indexing the bd
        % cell, returns a cell array data the same length as idx. In each
        % cell of data we have the following data per row:
        % [ length, normal, IDs ],
        % where length denotes the length of the respective edge on the
        % reference element, normal is the normals vector of the global
        % edge and IDs contains the indices of the edges' nodes.
        % Each cell of data corresponds to a requested sideset.
        %
            ss      = obj.bd(2,idx);
            Nss     = length(idx);
            meas    = cell(Nss);
            opp     = cell(Nss);
            normals = cell(Nss);
            faces   = cell(Nss);
            npere   = size(obj.el,2);

            switch obj.type
                case 'edge'
                    %% edges

                case 'tri'
                    %% triangles
                    % sideset coding for 2D triangle structures by the
                    % Abaqus inp file is as follows:
                    % E1: 1, 2
                    % E2: 2, 3
                    % E3: 3, 1
                    for i=1:Nss
                        switch npere
                            case 3
                                %% first order mesh
                                % 1 2 3
                                faces{i} = [obj.el(ss{1,i}{2,1},[1 2]);... % E1
                                            obj.el(ss{1,i}{2,2},[2 3]);... % E2
                                            obj.el(ss{1,i}{2,3},[3 1])];   % E3

                            case 6
                                %% second order mesh
                                % 1 2 3 e12 e23 e31
                                faces{i} = [obj.el(ss{1,i}{2,1},[1 2 4]);... % E1
                                            obj.el(ss{1,i}{2,2},[2 3 5]);... % E2
                                            obj.el(ss{1,i}{2,3},[3 1 6])];   % E3

                        end

                        opp{i}   = [obj.el(ss{1,i}{2,1},3); ... % E1
                                    obj.el(ss{1,i}{2,2},1); ... % E2
                                    obj.el(ss{1,i}{2,3},2)];    % E3

                        edges      = obj.co(:,:,faces{i}(:,2))-obj.co(:,:,faces{i}(:,1));
                        meas{i}    = sqrt(dot(edges,edges,1));
                        normals{i} = -edges.rot;
                        
                        % correct direction
                        tmp = ofem.matrixarray(ones(1,1,size(normals{i},3)));
                        tmp(1,1,dot(normals{i},obj.co(:,:,opp{i}),1)>=0) = -1;

                        normals{i} = normals{i}./(2*meas{i}.*tmp);
                    end
                            
                case 'quad'
                    %% quadrilateral
                    error('ofem:mesh:InvalidMesh',...
                          'Quadrilateral meshes not supported so far!');
                case 'tet'
                    %% tets
                    % sideset coding for 3D tetrahedral structures by
                    % the Abaqus inp file is as follows:
                    % S1: 1, 2, 3
                    % S2: 1, 2, 4
                    % S3: 2, 3, 4
                    % S4: 1, 3, 4
                    for i=1:Nss
                        switch npere
                            case 4
                                %% first order mesh
                                % 1 2 3 4
                                faces{i} = [obj.el(ss{1,i}{2,1},[1 2 3]);... % S1
                                            obj.el(ss{1,i}{2,2},[1 2 4]);... % S2
                                            obj.el(ss{1,i}{2,3},[2 3 4]);... % S3
                                            obj.el(ss{1,i}{2,4},[1 3 4])];   % S4

                            case 10
                                %% second order mesh
                                % 1 2 3 4 e12 e13 e14 e23 e24 e34
                                faces{i} = [obj.el(ss{1,i}{2,1},[1 2 3 5  8 6]);... % S1
                                            obj.el(ss{1,i}{2,2},[1 2 4 5  9 7]);... % S2
                                            obj.el(ss{1,i}{2,3},[2 3 4 8 10 9]);... % S3
                                            obj.el(ss{1,i}{2,4},[1 3 4 6 10 7])];   % S4

                        end

                        opp{i}   = [obj.el(ss{1,i}{2,1},4); ... % S1
                                    obj.el(ss{1,i}{2,2},3); ... % S2
                                    obj.el(ss{1,i}{2,3},1); ... % S3
                                    obj.el(ss{1,i}{2,4},2)];    % S4

                        e1         = obj.co(:,:,faces{i}(:,2))-obj.co(:,:,faces{i}(:,1));
                        e2         = obj.co(:,:,faces{i}(:,3))-obj.co(:,:,faces{i}(:,1));
                        normals{i} = cross(e1,e2,1);
                        meas{i}    = sqrt(dot(normals{i},normals{i},1))/2;

                        
                        % correct direction
                        tmp = ofem.matrixarray(ones(1,1,size(normals{i},3)));
                        tmp(1,1,dot(normals{i},obj.co(:,:,opp{i}),1)>=0) = -1;

                        normals{i} = normals{i}./(2*meas{i}.*tmp);
                    end
 
                case 'hex'
                    %% hexahedron
                    error('ofem:mesh:InvalidMesh',...
                          'Hexahedral meshes not supported so far!');
                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
        end

        function export_UCD(obj,folder_name,file_name,meta,varargin)
        %EXPORT_UCD exports solution to inp file.
        %
        % EXPORT_UCD(file_name,meta) exports a solution given by meta to a
        % inp file given by file_name. meta is explained below.
        %
        % EXPORT_UCD(file_name,meta1,meta2,...) exports several solutions.
        %
        % meta is a cell array numel(meta)==3, where meta{2} is the desired
        % function to be exported and must be known for every mesh node.
        % Note that ,e.g., for P2 elements the function must be known on
        % every additional node. meta{2} is not allowed to be sparse.
        % meta{1} is a string containing the name of the function to be
        % exported and meta{3} is a string containinig the name of the unit
        % associated to the function values.
        %
        % The exported inp file uses the UCD format, thus can be loaded,
        % e.g., in ParaView.
        %

            %% create output file
            if ~exist(folder_name,'dir')
                warning('ofem:mesh:FolderNotExistent',...
                        'The given folder is not existent, still I''ll create it.');
                mkdir(folder_name);
            end

            [pathstr,file_name,ext]=fileparts(file_name);
            if ~isempty(pathstr) % user requested current working directory
                warning('ofem:mesh:InvalidArgument',...
                        'file_name contains a path. I''ll ignore it');
            end
            if isempty(ext)
                file_name=strcat(file_name,'.inp');
            end
            file_name=fullfile(folder_name,file_name);


            %% basic computations
            Nf   = nargin-3         ; % number of functions
            Nc   = size(obj.el,1)   ; % number of cells
            Nn   = obj.Nco          ; % number of nodes
%             Nn  = size(obj.co,1);
            Nm   = size(obj.parts,2); % number of materials

            % 
            name = cell (Nf,1);
            unit = cell (Nf,1);
            Nfd  = zeros(Nf,1);

            if issparse(meta{2})
                warning('ofem:mesh:InvalidArgument',...
                        'u is not allowed to be sparse. I''ll convert it to full');
            end
            if numel(meta)~=3
                error('ofem:mesh:InvalidArgument',...
                      'meta is expected to be a cell array with numel(meta)==3');
            end
            name{1} = meta{1};
            sol     = full(meta{2});
            unit{1} = meta{3};
            Nfd (1) = size(meta{2},2);

            for i=1:nargin-4
                if issparse(varargin{i}{2})
                error('ofem:mesh:InvalidArgument',...
                      'u is not allowed to be sparse. I''ll convert it to full');
                end
                if numel(varargin{i})~=3
                    error('ofem:mesh:InvalidArgument',...
                          'meta is expected to be a cell array with numel(meta)==3');
                end

                name{i+1} = varargin{i}{1};
                Nfd (i+1) = size(varargin{i}{2},2);
                unit{i+1} = varargin{i}{3};

                sol(:,end+(1:Nfd(i+1))) = full(varargin{i}{2});
            end


            %% begin of export
            tic
            fprintf('Starting data export ... ');
            fileID = fopen(file_name,'w+');

            %% UCD format

            %% header => #nodes #cells #nodedata #celldata #classdata
            fprintf(fileID,'%d %d %d 0 0 \r\n',Nn, ...
                                               Nc, ...
                                               size(sol,2));

            %% nodes => nodeID <x_1> <x_2> <x_3>
            co_virt=zeros(Nn,3);
            co_virt(:,1:obj.dim)=squeeze(obj.co(:,1,1:Nn))';
            fprintf(fileID,'%d %e %e %e\r\n',[1:Nn;co_virt']);

            %% cells => cellID matID type <nodeID_1> ... <nodeID_N>
            currID=0;

            % As of now ParaView, or more precisely VTK, does not support
            % reading higher order elements although, at least, second
            % order elements are supported by UCD.
            switch obj.type
                case 'line'
                    formatStr = '%d %d line %d %d\r\n';
                    c2nidx = 1:2;
                case 'tri'
                    formatStr = '%d %d tri %d %d %d\r\n';
                    c2nidx = 1:3;
                case 'tet'
                    formatStr = '%d %d tet %d %d %d %d\r\n';
                    c2nidx = 1:4;
                case 'quad'
                    formatStr = '%d %d quad %d %d %d %d\r\n';
                    c2nidx = 1:4;
                case 'hex'
                    formatStr = '%d %d hex %d %d %d %d %d %d %d %d\r\n';
                    c2nidx = 1:8;
                otherwise
                    error('ofem:mesh:InvalidMesh',...
                          'Unsupported mesh found. Type is %s',obj.type);
            end
        
%             formatStr = ['%d %d ',obj.type,repmat(' %d',1,size(obj.el,2)),'\r\n'];
            
            for i=1:Nm
                partIDs  = obj.parts{3,i};
                NpartIDs = numel(partIDs);
                fprintf(fileID,formatStr,[currID+(1:NpartIDs); i*ones(1,NpartIDs); obj.el(partIDs,c2nidx)']);
                currID=currID+NpartIDs;
            end

            %% write functions
            % header => #func #colsfunc_1 ... #colsfunc_N
            fprintf(fileID,'%d ',[Nf; Nfd]);
            fprintf(fileID,'\r\n');

            % names, units
            for i=1:Nf
                fprintf(fileID,'%s, %s\r\n',name{i},unit{i});
            end

            % values => nodeID funcval_1 ... funcval_N
            formatStr = ['%d ',repmat(' %e',1,size(sol,2)),'\r\n'];
            fprintf(fileID,formatStr,[1:Nn;sol(1:Nn,:)']);

            %% close file print info and exit
            fclose(fileID);

            t=toc;
            fprintf('done. t=%f s\n',t);
        
        end

        %%
        function export_XDMF(obj,folder_name,file_name,meta,varargin)
        %EXPORT_XDMF exports solution to xdmf file.
        %
        % EXPORT_XDMF(file_name,meta) exports a solution given by meta to a
        % xml-file given by file_name. meta is explained below.
        %
        % EXPORT_XDMF(file_name,meta1,meta2,...) exports several solutions.
        %
        % meta is a cell array numel(meta)==3, where meta{2} is the desired
        % function to be exported and must be known for every mesh node.
        % Note that ,e.g., for P2 elements the function must be known on
        % every additional node. meta{2} is not allowed to be sparse.
        % meta{1} is a string containing the name of the function to be
        % exported and meta{3} is a string containinig the name of the unit
        % associated to the function values.
        %
        % The exported xml-file uses the XDMF format, thus can be loaded,
        % e.g., in ParaView.
        %

            error('ofem:mesh:NotImplemented',...
                  'XDMF file exporter is under development!');
        end

%         function info=printfinfo(obj)
%             edge = obj.el(:,[]
%             elem     = sortelem3(elem);
%             [~,edge] = dof3edge (elem);
%             [~,face] = dof3face (elem);
%             NE = size(edge,1); NF = size(face,1); NT = size(elem,1);
% 
%             % edge length statistics
%             e12 = node(edge(:,2),:)-node(edge(:,1),:);
%             edgeLength = sqrt(sum(e12.^2,2));
%             info.edge.num = NE;
%             info.edge.min = min(edgeLength);
%             info.edge.max = max(edgeLength);
%             info.edge.ave = sum(edgeLength)/NE;
% 
%             % areas
%             e12 = node(face(:,2),:)-node(face(:,1),:);
%             e13 = node(face(:,3),:)-node(face(:,1),:);
%             area = sqrt(sum(cross(e12,e13,2).^2,2))/2;
%             info.face.num = NF;
%             info.face.min = min(area);
%             info.face.max = max(area);
%             info.face.ave = sum(area)/NF;
% 
%             % volumes
%             e12 = node(elem(:,2),:)-node(elem(:,1),:);
%             e13 = node(elem(:,3),:)-node(elem(:,1),:);
%             e14 = node(elem(:,4),:)-node(elem(:,1),:);
%             volume = abs(dot(e12,cross(e13,e14,2),2))/6;
%             info.tet.num = NT;
%             info.tet.min = min(volume);
%             info.tet.max = max(volume);
%             info.tet.ave = sum(volume)/NT;
%         end

        %%
        function h=show(obj)
            switch obj.type
                case 'edge'
                    %% edge
                    warning('ofem:mesh:show:notImplemented',...
                            'Showing edge meshes not supported so far!');
 
                case 'tri'
                    %% triangle
                    h=trimesh(obj.el,obj.co(1,1,:),obj.co(2,1,:));
 
                case 'quad'
                    %% quadrilateral
                    warning('ofem:mesh:show:notImplemented',...
                            'Showing quadrilateral meshes not supported so far!');
 
                case 'tet'
                    %% tetrahedron
                    h=tetramesh(obj.el,double(permute(obj.co,[3,1,2])),'FaceAlpha',0.1);
 
                case 'hex'
                    %% hexahedron
                    warning('ofem:mesh:show:notImplemented',...
                            'Showing hexahedral meshes not supported so far!');
 
                otherwise
                    error('ofem:mesh:show:Unspecified',...
                          'Unspecified error found');
            end
        end
    end
end