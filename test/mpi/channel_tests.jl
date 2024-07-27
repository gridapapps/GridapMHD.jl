module ChannelTestsMPI

using GridapMHD: channel

function main(parts)
  channel(backend=:mpi,np=parts,sizes=(4,1,1),nc=(8,8,8,))
end

end # module
