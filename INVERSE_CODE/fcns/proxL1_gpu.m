function xhat_g = proxL1_gpu(y_g, lam)
xhat_g = y_g;
xhat_g = xhat_g + lam*((xhat_g<-lam)-(xhat_g>lam));
xhat_g(abs(xhat_g)<=lam)=0;

end