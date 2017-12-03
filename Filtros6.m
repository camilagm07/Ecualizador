function [k,N,D] = Filtros(SOS,G)

k = (G(1)*G(2)*G(3)*G(4))
N = conv([SOS(1) SOS(4) SOS(7)],conv([SOS(2) SOS(5) SOS(8)],[SOS(3) SOS(6) SOS(9)]))
D = conv([SOS(10) SOS(13) SOS(16)],conv([SOS(11) SOS(14) SOS(17)],[SOS(12) SOS(15) SOS(18)]))
end
