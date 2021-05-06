function xhat_g = Adam_proxL1_gpu(y_g, lam, phi)

xhat_g = y_g;
xhat_g = xhat_g - lam./phi.*(xhat_g>lam./phi) + lam./phi.*(xhat_g<-lam./phi);
xhat_g(abs(xhat_g)<=lam./phi)=0;

end