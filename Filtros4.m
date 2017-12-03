function [k,N,D] = Filtros4(SOS,G)
k = (G(1)*G(2)*G(3))
N = conv([SOS(1) SOS(3) SOS(5)],[SOS(2) SOS(4) SOS(6)])
D = conv([SOS(7) SOS(9) SOS(11)],[SOS(8) SOS(10) SOS(12)])
end