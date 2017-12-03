function [k,N,D] = Filtros8(SOS,G)

k = (G(1)*G(2)*G(3)*G(4)*G(5));
N = conv([SOS(1) SOS(5) SOS(9)],conv([SOS(2) SOS(6) SOS(10)],conv([SOS(3) SOS(7) SOS(11)],[SOS(4) SOS(8) SOS(12)])))
D = conv([SOS(13) SOS(17) SOS(21)],conv([SOS(14) SOS(18) SOS(22)],conv([SOS(15) SOS(19) SOS(23)],[SOS(16) SOS(20) SOS(24)])))
end
