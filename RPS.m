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

    % list of edges
    E = calculate_edges( vertices, N);
    
    % index of the next edge to be inserted
    edges_idx = 1;
    
    % Iterate through all the vertices to determine the visible vertices
    % from each vertex
    %for i=1:N
    for i=1:1
        
        % vertex v: start point
        v = vertices(i,:);

        % subset of vertex except the start point
        subset = [vertices(1:i-1,:) vertices(i+1:N,:)];

        % angle from the horizontal axis to the line segment vv_i sorted in
        % incresing order
        A = calculate_alpha( v, subset );

        % sorted list of edges that intersects the horizontal line
        % emanating from v
        S = intersects_line( v, vertices, E )
        
        % evaluate each vertex
        for j=1:N-1
                        
            vi = subset(A(i),:)
            
            if is_visible(v, vi)
                % add index of v and vi to the visibility graph
                edges(edges_idx) = [i, j+i];
            end
            
            %TODO: agregar los dos if
%             if 
%             end
%             
%             if
%             end
        end
        
    end
    
end

function [ E ] = calculate_edges( vertices, N )
%CALCULATE_EDGES Transforms the initial set of vertices into a data
%structure that represents the edges of the polygons.
%
%   E: Array of size M x 2, where M is the number of edges. Each row
%   represents one edge, having the start and the end vertices index of the
%   array of vertices.

    edge_idx = 1;

    for i=2:N-1
        
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
        alpha(i,1) = atan2( y, x ) ;
    end
    
    % converts angle to interval [0, 2*pi]
    alpha = mod(alpha,2*pi);
    
    % sort the elements of alpha
    [sorted_alpha, A] = sort(alpha(:),'ascend')
    
    % transpose A, so it has size 1 x N
    A = A';
    
end

function [ S ] = intersects_line( v, vertices, E )
%INTERSECTS_LINE Select the edges that intersects the horizontal half-line
%emanating from v
%
%   Detailed explanation goes here

    % number of edges
    N = size(E,1);

    % horizontal half-line emanating from v
    line_v = [v(1) v(2) vertices(end,1) v(2)];
    
    % index of the S array
    s_idx = 1;
    
    % determines the edges that intersects the line
    for i=1:N
    
        % define edge as a line in the form [x1 y1 x2 y2]
        line = [vertices( E(i,1), 1:2),  vertices( E(i,2), 1:2)];
        
        % determine wheter lines intersects or not and compute the distance
        % to the initial vertex
        [ intersect , dst] = is_intersected(line_v, line);
        
        % define whether the lines intersects or not
        if intersect
            % add edge index to the list
            S(s_idx) = i;
            % assign the distance
            E_dst(s_idx) = dst;
            s_idx = s_idx + 1;
        end
        
    end

    % sort the elements of E_dst and return the indexes
    [sorted_E_dst, sorted_idx] = sort(E_dst(:),'ascend');
    
    % sort S according to the sorted indexes
    S = S(sorted_idx');
    
end

function [ intersect, dst ] = is_intersected(line_1, line_2)
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
        
        % x-coordinate is on the lines
        if ( ( x > x_line_1(1) ) && ( x < x_line_1(2) ) && ...
             ( x > x_line_2(1) ) && ( x < x_line_2(2) ) )
            intersect = true;
            % euclidean distance
            dst = norm((intersect_pt(1:2) - line_1(1:2)),2);
        else
            intersect = false;
        end
        
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

function [ visible ] = is_visible( v, vi )
%IS_VISIBLE Summary of this function goes here
%   Detailed explanation goes here


end

% function [ intersect ] = is_intersected(line_1, line_2)
% %IS_INTERSECTED Determines whether two lines intersects or not
% %
% %   - line_1: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
% %     y1 is the initial point of the line and x2 y2 the end point
% %   - line_2: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
% %     y1 is the initial point of the line and x2 y2 the end point
% 
%     % transform line 1 to slope-intersect form
%     [m_1, b_1] = slope_intersect_form( line_1 );
%     
%     % transform line 2 to slope-intersect form
%     [m_2, b_2] = slope_intersect_form( line_2 );
%     
%     % determine if lines are parallel
%     if ( m_1 == m_2 )
%         intersect = false;
%     else
%         % x-coordinate of the intersection
%         x = ( b_2 - b_1 ) ./ ( m_1 - m_2 );
%         
%         % x-coordinate is on the lines
%         if ( ( x > line_1(1) ) && ( x < line_1(3) ) && ...
%              ( x > line_2(1) ) && ( x < line_2(3) ) )
%             intersect = true;
%         else
%             intersect = false;
%         end
%     end
%     
% end

% function [m, b] = slope_intersect_form( line )
% %SLOPE_INTERSECT_FORM Transforms a line having the starting and ending
% %points to the slope-insersect form
% %
% %   - line: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
% %     y1 is the initial point of the line and x2 y2 the end point
% 
%     m = ( line(4) - line(2) )./ ( line(3) - line(1) );
%     b = ( m .* -line(1) ) + line(2);
% 
% end


