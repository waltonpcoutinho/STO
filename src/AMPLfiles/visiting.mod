subject to
#"waypoint constraint"
waypoint_visit: (xBar - Y[N-1,1])^2 + (yBar - Y[N-1,2])^2 <= (Y[N-1,3] + rBar)^2;

#"photographing constraints"
levelH:     UbarH <= Y[N-1,3] <= barH;
levelG: -hatGamma <= Y[N-1,5] <= hatGamma;
levelM:    -hatMu <= U[N-1,2] <= hatMu;

