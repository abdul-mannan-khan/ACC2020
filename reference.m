function [pd,dot_pd,ddot_pd,d3dot_pd] = reference(t,varargin)
    pd = zeros(3,1);
    dot_pd = zeros(3,1);
    ddot_pd = zeros(3,1);
    d3dot_pd = zeros(3,1);
    if nargin > 1 && strcmpi(varargin{1},'circle')
        pd = [cos(t);sin(t);0];
        dot_pd = [-sin(t);cos(t);0];
        ddot_pd =[-cos(t);-sin(t);0];
        d3dot_pd =[sin(t);-cos(t);0];
    end
end