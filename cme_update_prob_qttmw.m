function v1 = cme_update_prob_qttmw(v,lb_old,l2size_old,lb_new,l2size_new)

    global no_species

    N = no_species;
    mode_old = 2.^l2size_old; % mode_old=[2^l1,...,2^lN]
    mode_new = 2.^l2size_new; % mode_new=[2^l'1,...,2^l'N]

    ub_old = lb_old + mode_old - 1;
    ub_new = lb_new + mode_new - 1;

    v_physic = tt_reshape(v,mode_old);
    v1_physic = change_prob(N,lb_old,ub_old,lb_new,ub_new,v_physic);
    v1 = tt_reshape(v1_physic,2*ones(sum(l2size_new),1));
end

function q = change_prob(N,lb_old,ub_old,lb_new,ub_new,p)
% INPUT :   N = no_species Global variable
%     0<=lb_old(s)<ub_old(s), 0<=lb_new(s)<ub_new(s),for s=1:N
%           p\in R^{I1*...*IN}, where Is=ub_old(s)-lb_old(s)+1
% OUTPUT : q\in R^{J1*...*JN}, where Js=ub_new(s)-lb_new(s)+1
%     so that ...
q = p; % Keep q.d=p.d and q.r=p.r
n_p = ub_old - lb_old + 1; % n_p(s)=Is
n_q = ub_new - lb_new + 1; % n_q(s)=Js, n_q will be q.n
% Prepare ps_q which will be q.ps
ps_q = p.ps; % ps_q has same size of p.ps and ps_q(1)=1
for i=2:N+1
      ps_q(i) = ps_q(i-1)+p.r(i-1)*n_q(i-1)*p.r(i);
end
% Prepare core_q which will be q.core
core_q = zeros(ps_q(N+1)-1,1);
for s=1:N
      core_s_p = p.core(p.ps(s):p.ps(s+1)-1);
      mat_s_p = reshape(core_s_p,[p.r(s),p.n(s),p.r(s+1)]);
      mat_s_q = zeros(p.r(s),n_q(s),p.r(s+1));
      move = lb_new(s) - lb_old(s);
      for i=1:n_q(s)
            k = i + move;
            if ((k >= 1) && (k <= n_p(s)))
                  mat_s_q(:,i,:) = mat_s_p(:,k,:);
            end
      end
      core_q(ps_q(s):ps_q(s+1)-1) = mat_s_q(:);
      core_q(ps_q(s):ps_q(s+1)-1);
end
q.n = n_q;
q.ps = ps_q;
q.core = core_q;
end
