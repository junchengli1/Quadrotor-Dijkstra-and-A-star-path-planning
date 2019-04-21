 function [path, num_expanded] = dijkstra(map, start, goal, astar)
% DIJKSTRA Find the shortest path from start to goal.
%   PATH = DIJKSTRA(map, start, goal) returns an mx3 matrix, where each row
%   consists of the (x, y, z) coordinates of a point on the path. The first
%   row is start and the last row is goal. If no path is found, PATH is a
%   0x3 matrix. Consecutive points in PATH should not be farther apart than
%   neighboring voxels in the map (e.g. if 5 consecutive points in PATH are
%   co-linear, don't simplify PATH by removing the 3 intermediate points).
%
%   PATH = DIJKSTRA(map, start, goal, astar) finds the path using euclidean
%   distance to goal as a heuristic if astar is true.
%
%   [PATH, NUM_EXPANDED] = DIJKSTRA(...) returns the path as well as
%   the number of nodes that were expanded while performing the search.
%   
% paramaters:
%   map     - the map object to plan in
%   start   - 1x3 vector of the starting coordinates [x,y,z]
%   goal:   - 1x3 vector of the goal coordinates [x,y,z]
%   astar   - boolean use astar or dijkstra
map = load_map('map1.txt', 0.1, 2, 0.2);
% % %%%for map 1%%%%%
start = [10  -4.9 2];
goal  = [6  18 2];
% start = [3  2 0];
% goal  = [5  2 4];
if nargin < 4
    astar = false;
end
start_sub=pos2sub(map,start);
goal_sub=pos2sub(map,goal);
start_ind=pos2ind(map,start);
goal_ind=pos2ind(map,goal);
num_expanded = 0;
path = [];
% Initialize distance array
[h w d]=size(map.occgrid);
%g[v]=inf
g=inf(h,w,d);
% priority queue of open nodes
Q=[];
Q=linspace(1,h*w*d,h*w*d);
%parent node
p=[];
H=zeros(h,w,d);
%distance from vs to vs is 0
g(start_ind)=0;
true;
while true
    if astar==true
        f=g+H;
    else
       f=g;
    end
[min_dist current_ind] = min(f(:));
 if ((current_ind == goal_ind) || isinf(min_dist))
        break
 end
[i,j,k]=ind2sub(size(g),current_ind);
%%%since we don't consider the diagonal, so we have 6 neighbours
neighbour =[i-1 j k;i j-1 k;i j k-1;i j k+1;i j+1 k; i+1 j k];
II=1;
Q(current_ind)=0;
%%%%%%%%check boundary and check collision
for ii=1:length(neighbour)
I=neighbour(:,1);
J=neighbour(:,2);
K=neighbour(:,3);
if (I(ii)<=h)&&(I(ii)>=1)&&(J(ii)<=w)&&(J(ii)>=1)&&(K(ii)<=d)&&(K(ii)>=1)...
   && map.occgrid(I(ii),J(ii),K(ii))==0
   neighbour_new(II,:)=neighbour(ii,:);
   II=II+1;   
end
end
neighbour=neighbour_new;
[I1,J1,K1]=ind2sub(size(map.occgrid),current_ind);
current_sub=[I1 J1 K1];
%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(neighbour)
neighbour_ind=sub2ind(size(map.occgrid),neighbour(ii,1),neighbour(ii,2),neighbour(ii,3));
    if Q(neighbour_ind)~=0
%         current_pos=ind2pos(map,current_ind);
        D=g(current_ind)+sum((map.res_xyz.*abs(neighbour(ii,:)-current_sub))');   
        if D<g(neighbour_ind)
            g(neighbour_ind)=D;
            p(neighbour_ind)=current_ind;
            H(neighbour_ind)=g(neighbour_ind)+norm(map.res_xyz.*abs(neighbour(ii,:)-goal_sub));
        end
    end
end
if current_ind~=goal_ind
        g(current_ind)=inf;
        num_expanded = num_expanded +1;
        end
end
if (isinf(g(goal_ind)))
path=[];
 else
    path=[goal_ind];
while (p(path(1))~=0)
    path= [p(path(1)),path];
end  
end
 path=ind2pos(map,path);
 path(1,:)=start;
 path(length(path),:)=goal;
end