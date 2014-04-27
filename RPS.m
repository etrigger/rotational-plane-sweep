function [ edges ] = RPS( vertices )
%RPS Program the Rotation Plane Sweep algorithm to compute the visibility 
%of the starting point with respect to all the vertices and goal point.
%
%   - vertices: array of size N x 3, where each row [x, y, n] represents
%     the x coordinate, y coordinate and object number. The first row
%     represents the START vertex and the last one the GOAL.
%   - edges: array of size M x 2, where each row [vi ve] represents the
%     initial vertex and the final vertex of the edge. M is the number of
%     edges, which is not pre-defined, given that it depends on the
%     problem.

    % number of vertices
    N = size(vertices,1);

    % intitialise the output vector
    edges = [];
    
    % index of the next edge to be inserted
    edges_idx = 1;
    
    % list of edges
    E = calculate_edges( vertices, N);
    
    % plot the initial environment
    figure;
    plot_environmnet( E, vertices );
    
    % test whether the environment has no obstacles
    if( isempty(E) ); edges = [1,2]; plot_result( edges, vertices ); return; end;
    
    % Iterate through all the vertices to determine the visible vertices
    % from each vertex
    for i=1:N
        
        % vertex v: start point
        v = vertices(i,:);

        % subset of vertex except the start point
        subset = [vertices(1:i-1,:); vertices(i+1:N,:)];

        % angle from the horizontal axis to the line segment vv_i sorted in
        % incresing order
        A = calculate_alpha( v, subset );

        % sorted list of edges that intersects the horizontal line
        % emanating from v
        [S, E_dst] = intersects_line( v, vertices, E )
        
        % evaluate each vertex
        for j=1:N-1
            
            % determine the index of the vertex in the initial array
            if (A(j)<i); vertex_nr = A(j); else vertex_nr = A(j) + 1; end;
            
            % vertex whose visibility will be tested
            vi = subset(A(j),:);
            
            % add the edge to the visible list, in case is visible
            [edges, edges_idx] = add_edge(S, v, vi, E, vertices, E_dst, ...
                                          vertex_nr, edges_idx, edges, i);
                                 
            % determine the edges indexes where vi is the start edge
            start_edge = find( E(:,1) == vertex_nr );
            
            % determine the edges indexes where vi is the end edge
            end_edge = find( E(:,2) == vertex_nr );
            
            % find the edges that should be either deleted or inserted
            [insert_edges, delete_edges] = find_edges(v,vi,start_edge,end_edge,E,vertices);
            
            % if vi is in the begining of an edge that is not in S
            if ~isempty( insert_edges)
                % insert the edge in S
                [S, E_dst] = insert_edge(v, vi, insert_edges, E_dst, S, vertices);
            end
            
            % if vi is in the end of an edge in S
            if ~isempty( delete_edges)
                % delete the edge from S
                [S, E_dst] = delete_edge(delete_edges, E, E_dst, S);
            end
        end   
    end
    
    % delete internal edges and add polygon lines
    edges = clear_edges(edges, vertices, E);
    
    % plots the resulting edges
    plot_result( edges, vertices );
    % paint in blue the polygons
    plot_environmnet( E, vertices );
    
end

function [edges, edges_idx] = add_edge(S, v, vi, E, vertices, E_dst, vertex_nr, edges_idx, edges, i)
%ADD_EDGE Add visible edges to the vector edges and update the value of the
%index

    if ( ~isempty( S) )
    % test whether the vertex is visible
        if is_visible(v, vi, S, E, vertices, E_dst)
            % add indexes [v vi] to the visibility graph
            edges(edges_idx,:) = [i, vertex_nr];                    
            edges_idx = edges_idx + 1;
        end
    else
            % if S is empty, add index to the visibility graph
            edges(edges_idx,:) = [i, vertex_nr];
            edges_idx = edges_idx + 1;
    end
                                 
end

function [ insert_edges, delete_edges ] = find_edges( v, vi, start_idx, end_idx, E, vertices)
%FIND_EDGES Determines the initial and the end edges

    insert_edges = [];
    delete_edges = [];

    proj_start = 0;
    proj_end = 0;
    
    line_v_vi_hmg =  homogeneous_coordinates( [v(1) v(2) vi(1) vi(2)] );
    line_v_vi_hmg = line_v_vi_hmg ./ sqrt( line_v_vi_hmg(1).^2 + line_v_vi_hmg(2).^2);
    
    vertex_start = [vertices( E(start_idx,2), 1:2), 1];
    vertex_end = [vertices( E(end_idx,1), 1:2), 1];
   
    if( ~isempty(start_idx) )
        proj_start = dot(line_v_vi_hmg, vertex_start);
    end
    
    if( ~isempty(end_idx) )
        proj_end = dot(line_v_vi_hmg, vertex_end)
    end

    insert_idx = 1;
    delete_idx = 1;
    
    if( proj_start > 0  )
        insert_edges( insert_idx ) = start_idx;   
        insert_idx = insert_idx + 1;
    end
    if( proj_start < 0 ) 
        delete_edges( delete_idx ) = start_idx;   
        delete_idx = delete_idx + 1;
    end
    
    if( proj_end > 0  )
        insert_edges( insert_idx ) = end_idx;   
        insert_idx = insert_idx + 1;
    end
    if( proj_end < 0 ) 
        delete_edges( delete_idx ) = end_idx;   
        delete_idx = delete_idx + 1;
    end

end

function [ E ] = calculate_edges( vertices, N )
%CALCULATE_EDGES Transforms the initial set of vertices into a data
%structure that represents the edges of the polygons.
%
%   E: Array of size M x 2, where M is the number of edges. Each row
%   represents one edge, having the start and the end vertices index of the
%   array of vertices.

    % initialization of the variables
    E =[];
    edge_idx = 1;

    for i=2:N-1
        
        % number of the current polygon
        object_nr = vertices(i,3);
        
        % check previous vertex to find the initial vertex of the object
        if ( vertices(i-1,3) ~= object_nr )
            % initial vertex index of this object
            init_vertex_idx = i;
        end
        
        % check next vertex to find an edge
        if ( vertices(i+1,3) == object_nr )
            % the next vertex belongs to the same object, so there is an
            % edge
            E(edge_idx,:) = [i, i+1];
        else
            % the current vertex is the last one of the object, so there is
            % an edge between this vertex and the initial vertex of the
            % object
            E(edge_idx,:) = [i, init_vertex_idx];
        end
        
        edge_idx = edge_idx + 1;
        
    end
    
end


function [ A ] = calculate_alpha( v, vertices )
%CALCULATE_ALPHA Calculates the angle from the horizontal axis, to the line
%segment defined by each vertex of the vertices array and v. 
%
%   - vertices: array of size N x 3, where each row [x, y, n] represents
%     the x coordinate, y coordinate and object number. 
%   - v: start vertex, that is not contained on the vertices array.

    % number of vertices 
    N = size(vertices,1);
    
    % calculates the angle from the horizontal axis
    for i=1:N
        x = vertices(i,1) - v(1);
        y = vertices(i,2) - v(2);
        % computes the angle in the interval [-pi,pi]
        alpha(i,1) = atan2( y, x );
    end
    
    % converts angle to interval [0, 2*pi]
    alpha = mod(alpha,2*pi);
    
    % sort the elements of alpha
    [sorted_alpha, A] = sort(alpha(:),'ascend');
    
    % transpose A, so it has size 1 x N
    A = A';
    
end

function [ S, sorted_E_dst ] = intersects_line( v, vertices, E )
%INTERSECTS_LINE Select the edges that intersects the horizontal half-line
%emanating from v

    % initialisation of the vectors
    E_dst = [];
    S = [];

    % number of edges
    N = size(E,1);

    % maximum x-coordinate
    [~, max_idx] = max(vertices(:,1));
    
    % horizontal half-line emanating from v
    line_v = [v(1) v(2) vertices(max_idx,1) v(2)];
    
    % index of the S array
    s_idx = 1;
    
    % determines the edges that intersects the line
    for i=1:N
    
        % define edge as a line in the form [x1 y1 x2 y2]
        line = [vertices( E(i,1), 1:2),  vertices( E(i,2), 1:2)];
        
        % determine wheter lines intersects or not and compute the distance
        % to the initial vertex
        [ intersect , dst] = is_intersected(line_v, line, E, vertices);
        
        % define whether the lines intersects or not
        if intersect
            
            % determine if the lines are different segments of the same
            % line. If this is the case, the edge is not occluding any
            % line, so it is not included on the S list
            if abs( homogeneous_coordinates(line_v) - homogeneous_coordinates(line) ) > 0.00001
                % add edge index to the list
                S(s_idx) = i;
                % assign the distance
                E_dst(s_idx) = dst;
                s_idx = s_idx + 1;
            end
      
        end
        
    end

    % sort the elements of E_dst and return the indexes
    [sorted_E_dst, sorted_idx] = sort(E_dst(:),'ascend');
    sorted_E_dst = sorted_E_dst';

    % sort S according to the sorted indexes
    S = S(sorted_idx');
  
end

function [ intersect, dst ] = is_intersected(line_1, line_2, E, vertices)
%IS_INTERSECTED Determines whether two lines intersects or not
%
%   - line_1: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
%     y1 is the initial point of the line and x2 y2 the end point
%   - line_2: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
%     y1 is the initial point of the line and x2 y2 the end point

    % initialise the value of the distance
    dst = 0;

    % transform line 1 to homogeneous coordinates
    line_hmg_1 = homogeneous_coordinates( line_1 );
    
    % transform line 2 to slope-intersect form
    line_hmg_2 = homogeneous_coordinates( line_2 );
    
    % intersection point
    intersect_pt = cross(line_hmg_1,line_hmg_2);
    
    if intersect_pt(3) == 0
        intersect = false;
    else
        % normalize the vector, so the third component is one
        intersect_pt = intersect_pt ./ intersect_pt(3);
        
        % x-coordinate of the intersection
        x = intersect_pt(1);
        
        % sort the x values of each line
        x_line_1 = sort([line_1(1), line_1(3)],'ascend');
        x_line_2 = sort([line_2(1), line_2(3)],'ascend');
        
        %check if the intersection is on an edge. In this case, there is no
        %occlusion
        if ( is_an_edge( intersect_pt, E, vertices, line_1) || ...
             is_an_edge( intersect_pt, E, vertices, line_2)  )
            intersect = false;
        else
            % x-coordinate is on the lines
            if ( ( x >= x_line_1(1) ) && ( x <= x_line_1(2) ) && ...
                 ( x >= x_line_2(1) ) && ( x <= x_line_2(2) ) )
                intersect = true;
                % euclidean distance
                dst = norm((intersect_pt(1:2) - line_1(1:2)),2);
            else
                intersect = false;
            end
        end
        
    end
    
end

function [ edge ] = is_an_edge( intersect_pt, E, vertices, line )
%IS_AN_EDGE Determines whether the point defined by intersect_pt lies in
%one of the vertex of and edge.
%
% The difference between both coordinates is computed and evaluated agains
% a threshold.

     diff_x_line_1 = abs( intersect_pt(1) - line(1) );
     diff_y_line_1 = abs( intersect_pt(2) - line(2) );
     diff_x_line_2 = abs( intersect_pt(1) - line(3) );
     diff_y_line_2 = abs( intersect_pt(2) - line(4) );
     
     if( ( ( diff_x_line_1 < 0.0001 ) && ( diff_y_line_1 < 0.0001 ) ) || ...
         ( ( diff_x_line_2 < 0.0001 ) && ( diff_y_line_2 < 0.0001 ) ) )
        edge = true;
     else
         edge = false;
     end

end

function [hmg] = homogeneous_coordinates( line )
%HOMOGENEOUS_COORDINATES Transforms a line having the starting and ending
%points to homogeneous coordinates
%
%   - line: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
%     y1 is the initial point of the line and x2 y2 the end point

    a = line(2) - line(4);
    b = line(3) - line(1);
    c = ( line(1) * line(4) ) - ( line(3) * line(2) );

    hmg = [a, b, c];
end

function [ visible ] = is_visible(v, vi, S, E, vertices, E_dst)
%IS_VISIBLE Determines whether a vertex is visible or not

    % distance from v to vi
    dst = norm((v(1:2) - vi(1:2)),2);
    
    % number of edges in S
    N = size(S,2);
    
    % index of S
    S_idx = 1;
    
    % initialise the visible variable
    visible = true;
    
    % defines the line that goes from v to the evaluated vertex, vi
    line_v_vi = [v(1) v(2) vi(1) vi(2)];
    
    % determine if the line v_vi intersects with any edge that is closer to
    % the candidate vertex
%     while ( ( S_idx <= N ) && ( dst >= E_dst(S_idx) ) )
    while ( ( S_idx <= N ) )
        
        E_idx = S(S_idx)
        
        line_e = [vertices( E(E_idx,1), 1:2),  vertices( E(E_idx,2), 1:2)];
        
        [ intersect , ~] = is_intersected(line_v_vi, line_e, E, vertices);
        
        if( intersect )
            visible = false;
            break;
        end
        
        S_idx = S_idx + 1;
        
    end

end

function [ S, E_dst ] = insert_edge(v, vi, edges_insrt, E_dst, S, vertices)
%INSERT_EDGES Summary of this function goes here
%   Detailed explanation goes here

    % number of edges to be inserted
    N = size(edges_insrt,2);
    
    for i=1:N
        
        % the edge is not in S
        if isempty( find(S(:) == edges_insrt(i) ) )
        
            % TODO: correct this distance
        
            % distance from the origin to the vertex of the edge
            dst = norm((v(1:2) - vi(1:2)),2);

            % TODO: fix this =, because is an special case
            smaller_idx = E_dst(:) <= dst;
            bigger_idx = E_dst(:) > dst;
            
            % insert edge into S
            S = [S(smaller_idx'), edges_insrt(i), S(bigger_idx') ];
            
            % insert distance into E_dst
            E_dst = [E_dst(smaller_idx'), dst, E_dst(bigger_idx')];
        end
            
    end

end

function [ S_out, E_dst_out ] = delete_edge(edges_dlt, E, E_dst, S);
%INSERT_EDGES Summary of this function goes here
%   Detailed explanation goes here

    % number of edges to be deleted
    N = size(edges_dlt,2);
    
    S_size = size(S,2);
    
    S_out = S;
    E_dst_out = E_dst;
    
    for i=1:N
        
        % the edge is not in S
        if ~( isempty( find(S_out(:) == edges_dlt(i)) ) )
            
            idx_S = find( S_out(:) == edges_dlt(i) );
            
            S_out(idx_S) = [];
            E_dst_out(idx_S) = [];
        
%             % insert edge into S
%             S_out = [S(1:idx_S - 1), S(idx_S + 1:S_size)];
% 
%             % insert distance into E_dst
%             E_dst_out = [E_dst(1:idx_S - 1), E_dst(idx_S + 1:S_size)];

        end
        
    end

end

function [edges] = clear_edges(edges, vertices, E)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    edge_idx = 1;
    
    % delete edges that belongs to the same polygon
    for i=1:size(edges,1)
        if ( vertices( edges(edge_idx,1),3 ) == vertices( edges(edge_idx,2),3 ) )
            edges(edge_idx,:) = [];
        else
            edge_idx = edge_idx + 1;
        end
    end
    
    edges = [edges; E];
    edges = sortrows(edges);

end

function plot_environmnet( E, vertices )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    % start 
    plot( vertices(1,1), vertices(1,2), 'ro', 'MarkerFaceColor',[1,0,0], ...
          'MarkerSize', 8 );
    hold on;
    % goal
    plot( vertices(end,1), vertices(end,2), 'go', 'MarkerFaceColor',[0,1,0], ...
          'MarkerSize', 8 );
    hold on;
    
    % Plot the edges
    for i=1:size(E,1)
        plot([vertices( E(i,1), 1); vertices( E(i,2), 1)],...
             [vertices( E(i,1), 2); vertices( E(i,2), 2)],...
             'bo', 'MarkerFaceColor',[0,0,1]);
        hold on;
        plot([vertices( E(i,1), 1); vertices( E(i,2), 1)],...
             [vertices( E(i,1), 2); vertices( E(i,2), 2)]);
        hold on;
    end
    
end

function plot_result( edges, vertices )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    for n=1:size(edges,1)
        plot([vertices( edges(n,1), 1); vertices( edges(n,2), 1)],...
             [vertices( edges(n,1), 2); vertices( edges(n,2), 2)],...
             '-r');
        hold on;
    end
    
end

function [ vertices ] = environment( )
%ENVIRONMENT Construct an environment based on the input provided by the
%user. The number polygons constructed must have four sides, and any number
%of them could be built.
%
%   - vertices: array [x y n], where x and y represents the coordinate of
%   the point and n is the number of the polygon to which it belongs.

    vertices = zeros(10,3);

    [x,y] = ginput(10);

    vertices(1,:) = [x(1) y(1) 0];
    vertices(2:5,:) = [x(2:5) y(2:5) ones(4,1)];
    vertices(6:9,:) = [x(6:9) y(6:9) (ones(4,1).*2)];
    vertices(10,:) = [x(10) y(10) 3];

    plot( x(1), y(1), 'bo' )
    hold on;
    plot([x(2:5); x(2)],[y(2:5); y(2)]);
    hold on;
    plot([x(6:9); x(6)],[y(6:9); y(6)]);
    hold on;
    plot(x(10),y(10), 'ro');

end
