function [mid_point,len,vector] = segment_calc(P1,P2)
%SEGMENT_CALC calculates the middle point of a segment given the two
%extreme points of the segment. it also calculates the length of the
%segment and the vector pointing from one of the extreme points of the
%segment to the other.
%   P1 and P2 are the coordinates of the two extreme points of the segment.
vector = P2-P1;
len = norm(vector);
mid_point = P1 + (vector/2);
end

