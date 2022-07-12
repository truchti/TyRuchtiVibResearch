kv1 = open_loop_knot_vector(20, 0, 4)
kv2 = custom_boundary_open_knot_vectors(20, 0, 4, 3,5)

function knotVector = open_loop_knot_vector(numPts, minVal, maxVal)
npc = numPts+ 6;
np2 = numPts+2;
knotVector = zeros(1,npc); % total number of knots is number of control points plus the order in this case 6
for i = 2:npc
    if i>6 && i<np2
        knotVector(i) = knotVector(i-1) + 1;
    else
        knotVector(i) = knotVector(i-1);
    end
end
if minVal < 0
    knotVector = (knotVector - max(knotVector))./max(knotVector)*abs(minVal);
else
    knotVector = knotVector./max(knotVector)*maxVal;
end
end
function knotVector = custom_boundary_open_knot_vectors(numPts, minVal, maxVal, startRepeat, endRepeat)
npc = numPts + 6;
np2 = numPts + 2;
knotVector = zeros(1,npc);
for i = 2:npc
    if i > startRepeat && i < (np2 + (6-endRepeat))
        knotVector(i) = knotVector(i-1) + 1;
    else
        knotVector(i) = knotVector(i-1);
    end
end
if minVal < 0
    knotVector = (knotVector - max(knotVector))./max(knotVector)*abs(minVal);
else
    knotVector = knotVector./max(knotVector)*maxVal;
end
end