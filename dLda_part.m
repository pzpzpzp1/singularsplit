function [val, grad, extra] = dLda_part(dAda, fun, ls, ehat, x0, extraIn)
    if exist('extraIn','var')
        u1 = extraIn.u1;
        v1 = extraIn.v1;
        u3 = extraIn.u3;
        v3 = extraIn.v3;
        Eic = extraIn.Eic;
    else
        u1 = [0;0];
        v1 = [0;0];
        u3 = [0;0];
        v3 = [0;0];
        Eic = [0 0];
    end


    uv45 = x0;
    u4 = uv45(1:2);
    u5 = uv45(3:4);
    v4 = uv45(5:6);
    v5 = uv45(7:8);
    l1u = ls(1);
    l1v = ls(2);
    l2u = ls(3);
    l2v = ls(4);
    [E4,dEduv4] = fun([u4;v4]);
    [E5,dEduv5] = fun([u5;v5]);

    val1 = dot(dAda, [Eic E4 E5]);
    val2 = -(l1u*dot(u1-u4,ehat) + l1v*dot(v1-v4,ehat) + l2u*dot(u3-u5,ehat) + l2v*dot(v3-v5,ehat));
    val = val1+val2;
    
    z = [0 0];
    duv4 = dEduv4*dAda(3) + (l1u*[ehat; z'] + l1v*[z'; ehat]);
    duv5 = dEduv5*dAda(4) + (l2u*[ehat; z'] + l2v*[z'; ehat]);
    duv45 = reshape([reshape(duv4,2,2); reshape(duv5,2,2); ],[],1);
    grad = duv45;
    
    extra.val1 = val1;
    extra.val2 = val2;
    extra.dEduv4 = dEduv4;
    extra.dEduv5 = dEduv5;
end