
% state the current weather as a list of facts:

outlook(overcast).   % can be either sunny, overcast or rain
humidity(high).      % humidity can be either normal or high
wind(strong).        % wind can be strong or weak

%---- solution 1:
playtennis :- outlook(overcast).
playtennis :- outlook(sunny), humidity(normal).
playtennis :- outlook(rain), wind(weak).

%---- solution 2:
playtennis :- outlook(overcast) ;
              outlook(sunny), humidity(normal) ;
              outlook(rain), wind(weak).

% now change the outlook, humidity and wind to the allowed values and
% ask prolog if you should play tennis or not:

?- playtennis.
