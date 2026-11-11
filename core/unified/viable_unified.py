
from fractions import Fraction as Q
from .rational_polytope import Poly,box2,reach,residual

def data(full_disturbance=True,eta=None):
 Ai=((Q(1),Q(-1)),(Q(1),Q(-1)))
 Aj=((Q(41,40),Q(-21,20)),(Q(1),Q(-1)))
 B=((Q(1,10),Q(0)),(Q(0),Q(1,10)))
 X=box2(-1,1,-1,1);X0=Poly([(Q(1),Q(0))]);U=box2(-1,1,-1,1)
 eta=Q(0) if not full_disturbance else (Q(1,10000) if eta is None else Q(eta))
 Wi=box2(-eta,eta,-eta,eta)
 Wj=box2(-Q(1,40)-eta,Q(11,400)+Q(41,40)*eta,-eta,eta)
 C=U.affine(B).minkowski(U.affine(((Q(0),Q(-1,10)),(Q(1,10),Q(-1,10)))))
 AB=((Q(1,10),Q(-1,10)),(Q(1,10),Q(-1,10)))
 C=U.affine(B).minkowski(U.affine(AB)).minkowski(Wi).minkowski(Wi.affine(Ai))
 return Ai,Aj,B,X,X0,U,Wi,Wj,C

def pipeline(H=2,full_disturbance=True,eta=None):
 Ai,Aj,B,X,X0,U,Wi,Wj,C=data(full_disturbance,eta=eta)
 safe=X0
 stages=[]
 for t in range(H):
  res=residual(safe,Ai,Aj,B,B,U,Wi)
  accepted=Wj.contains_poly(res)
  next_all=reach(safe,Ai,B,U,Wi)
  next_safe=next_all.clip_box()
  stages.append({'t':t,'safe_predecessors':safe.encode(),'source_successors_all':next_all.encode(),
                 'first_exit_residual':res.encode(),'representative_facet_pass':accepted,
                 'safe_successor_projection':next_safe.encode()})
  if not accepted:
   return {'status':'REJECT','finite_stages':stages,'tail':None}
  safe=next_safe
 tail_image=reach(C,Ai,B,U,Wi)
 tail_res=residual(C,Ai,Aj,B,B,U,Wi)
 capture=C.contains_poly(safe)
 safe_guard=X.contains_poly(C)
 closure=C.contains_poly(tail_image.clip_box())
 closure_full=C.contains_poly(tail_image)
 resid=Wj.contains_poly(tail_res)
 status='CERTIFIED' if capture and safe_guard and closure and resid and not safe.empty else 'UNKNOWN'
 return {'status':status,'horizon':H,'finite_stages':stages,
         'terminal_safe_reachable':safe.encode(),'terminal_safe_nonempty':not safe.empty,
         'tail':{'candidate_C':C.encode(),'capture':capture,'within_safety':safe_guard,
                 'safe_successor_closed':closure,'all_action_invariant':closure_full,
                 'residual_included':resid,'image':tail_image.encode(),
                 'tail_residual':tail_res.encode()},
         'model_quotient':{'original_models':2,'retained_models':1,'fixed_representative':'j'} }
