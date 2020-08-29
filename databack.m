function [] = databack(hlines)
%% Set the maximum number of lines to 12
number = length(hlines);
if number == 1
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
elseif number == 2
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
elseif number == 3
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
elseif number == 4
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
elseif number == 5
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(5),@(p) eval('uiresume'));
elseif number == 6
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(5),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(6),@(p) eval('uiresume'));
elseif number == 7
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(5),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(6),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(7),@(p) eval('uiresume'));
elseif number == 8
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(5),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(6),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(7),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(8),@(p) eval('uiresume'));
elseif number == 9
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(5),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(6),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(7),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(8),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(9),@(p) eval('uiresume'));
elseif number == 10
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(5),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(6),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(7),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(8),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(9),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(10),@(p) eval('uiresume'));
elseif number == 11
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(5),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(6),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(7),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(8),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(9),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(10),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(11),@(p) eval('uiresume'));
elseif number == 12
    addNewPositionCallback(hlines(1),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(2),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(3),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(4),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(5),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(6),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(7),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(8),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(9),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(10),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(11),@(p) eval('uiresume'));
    addNewPositionCallback(hlines(12),@(p) eval('uiresume'));
else
    error('The number of lines exceeds the limit.')
end
