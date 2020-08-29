function [] = databack(hlines)
%% Set the maximum number of lines to 12
number = length(hlines);
if number == 1
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
elseif number == 2
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
elseif number == 3
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
elseif number == 4
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
elseif number == 5
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(5),@(p) evalin('uiresume'));
elseif number == 6
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(5),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(6),@(p) evalin('uiresume'));
elseif number == 7
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(5),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(6),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(7),@(p) evalin('uiresume'));
elseif number == 8
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evailin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(5),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(6),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(7),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(8),@(p) evalin('uiresume'));
elseif number == 9
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(5),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(6),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(7),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(8),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(9),@(p) evalin('uiresume'));
elseif number == 10
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(5),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(6),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(7),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(8),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(9),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(10),@(p) evalin('uiresume'));
elseif number == 11
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(5),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(6),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(7),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(8),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(9),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(10),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(11),@(p) evalin('uiresume'));
elseif number == 12
    addNewPositionCallback(hlines(1),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(2),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(3),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(4),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(5),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(6),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(7),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(8),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(9),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(10),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(11),@(p) evalin('uiresume'));
    addNewPositionCallback(hlines(12),@(p) evalin('uiresume'));
else
    error('The number of lines exceeds the limit.')
end