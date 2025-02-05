#ifndef MPI_SERIALIZABLE_COMMANDS_HPP
#define MPI_SERIALIZABLE_COMMANDS_HPP

#ifdef MADVORO_WITH_MPI

#include <functional>
#include <mpi.h>
#include "Serializer.hpp"
#include "exception/MadVoroException.hpp"

#define MPI_EXCHANGE_TAG 5

namespace MadVoro
{
	using rank_t = int;

	namespace MPI
	{
		template<typename T, template<typename...> class Container, typename... Ts>
		std::vector<std::vector<T>> MPI_Exchange_all_to_all(const std::vector<Container<T, Ts...>> &data, const MPI_Comm &comm)
		{
			rank_t size;
			Serializer send;
			MPI_Comm_size(comm, &size);

			std::vector<int> sendDisplacements(size, 0), recvDisplacements(size, 0);
			std::vector<int> sendCounts(size, 0), recvCounts(size, 0);

			for(rank_t _rank = 0; _rank < size; _rank++)
			{
				sendCounts[_rank] = static_cast<int>(send.insert_all(data[_rank]));
				if(_rank > 0)
				{
					sendDisplacements[_rank] = sendDisplacements[_rank - 1] + sendCounts[_rank - 1];
				}
			}

			int totalSize = 0;
			MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, comm);

			for(rank_t _rank = 0; _rank < size; _rank++)
			{
				totalSize += recvCounts[_rank];
				if(_rank > 0)
				if(_rank > 0)
				{
					recvDisplacements[_rank] = recvDisplacements[_rank - 1] + recvCounts[_rank - 1];
				}
			}

			Serializer recv;
			recv.resize(totalSize);

			MPI_Alltoallv(send.getData(), sendCounts.data(), sendDisplacements.data(), MPI_BYTE, recv.getData(), recvCounts.data(), recvDisplacements.data(), MPI_BYTE, comm);

			std::vector<std::vector<T>> result(size);
			for(rank_t _rank = 0; _rank < size; _rank++)
			{
				size_t bytesRead = recv.extract(result[_rank], recvDisplacements[_rank], recvCounts[_rank]);
				assert(bytesRead != static_cast<size_t>(recvCounts[_rank]));
			}

			return result;
		}

		template<typename T>
		std::vector<std::vector<T>> MPI_Exchange_by_ownership_by_ranks(const std::vector<T> &data, const std::function<rank_t(const T&)> &ownership, const MPI_Comm &comm)
		{
			rank_t rank, size;
			MPI_Comm_rank(comm, &rank);
			MPI_Comm_size(comm, &size);

			std::vector<std::vector<T>> dataByRanks(size);
			for(const T &value : data)
			{
				rank_t _rank = ownership(value);
				if(_rank < 0 or _rank >= size)
				{
					MadVoro::Exception::MadVoroException eo("MPI_Exchange_by_ownership_by_ranks: ownership function returned invalid rank");
					eo.addEntry("Rank", _rank);
					eo.addEntry("Size", size);
					throw eo;
				}
				dataByRanks[_rank].push_back(value);
			}

			return MPI_Exchange_all_to_all(dataByRanks, comm);
		}

		template<typename T>
		std::vector<T> MPI_Exchange_by_ownership(const std::vector<T> &data, const std::function<rank_t(const T&)> &ownership, const MPI_Comm &comm)
		{
			std::vector<std::vector<T>> exchangedData = MPI_Exchange_by_ownership_by_ranks(data, ownership, comm);
			std::vector<T> result;
			for(const std::vector<T> &values : exchangedData)
			{
				result.insert(result.end(), values.cbegin(), values.cend());
			}
			return result;
		}

		template<typename T, template<typename...> class Container, typename... Ts>
		std::vector<std::vector<T>> MPI_All_cast_by_ranks(const Container<T, Ts...> &data, const MPI_Comm &comm = MPI_COMM_WORLD)
		{
			rank_t size;
			MPI_Comm_size(comm, &size);

			// first know how much data is being sent from each one
			Serializer send;
			int count = static_cast<int>(send.insert_all(data));
			std::vector<int> recvCounts(size, 0);
			MPI_Allgather(&count, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, comm);

			std::vector<int> recvDisplacements(size, 0);
			size_t totalToReceive = 0;
			for(rank_t _rank = 0; _rank < size; _rank++)
			{
				totalToReceive += static_cast<size_t>(recvCounts[_rank]);
				if(_rank > 0)
				{
					recvDisplacements[_rank] = recvDisplacements[_rank - 1] + recvCounts[_rank - 1];
				}
			}

			std::vector<int> sendDisplacements(size, 0);
			std::vector<int> sendCounts(size, count);
			Serializer recv;
			recv.resize(totalToReceive);
			MPI_Alltoallv(send.getData(), sendCounts.data(), sendDisplacements.data(), MPI_BYTE,
							recv.getData(), recvCounts.data(), recvDisplacements.data(), MPI_BYTE, comm);

			std::vector<std::vector<T>> resultByRanks(size);
			for(rank_t _rank = 0; _rank < size; _rank++)
			{
				size_t readCount = recv.extract(resultByRanks[_rank], recvDisplacements[_rank], recvCounts[_rank]);
				assert(readCount == recvCounts[_rank]);
			}
			return resultByRanks;
		}

		template<typename T, template<typename...> class Container, typename... Ts>
		std::vector<T> MPI_All_cast(const Container<T, Ts...> &data, const MPI_Comm &comm)
		{
			std::vector<std::vector<T>> resultByRanks = MPI_All_cast_by_ranks(data, comm);
			std::vector<T> result;
			for(const std::vector<T> &values : resultByRanks)
			{
				result.insert(result.end(), values.cbegin(), values.cend());
			}
			return result;
		}

		template<typename T>
		T MPI_Bcast_serializable(const T &data, rank_t owner, const MPI_Comm &comm = MPI_COMM_WORLD)
		{
			Serializer send;
			size_t sizeSent = send.insert(data);
			MPI_Bcast(&sizeSent, 1, MPI_UNSIGNED_LONG_LONG, owner, comm);
			
			Serializer recv;
			recv.resize(sizeSent);
			MPI_Bcast(recv.getData(), sizeSent, MPI_BYTE, owner, comm);

			T value;
			size_t sizeRead = recv.extract(value, 0);
			assert(sizeRead == sizeSent);
			return value;
		}

		template<typename T, template<typename...> class Container, typename... Ts>
		std::vector<T> MPI_Spread(const Container<T, Ts...> &data, rank_t root, const MPI_Comm &comm)
		{
			rank_t rank, size;
			MPI_Comm_rank(comm, &rank);
			MPI_Comm_size(comm, &size);

			if(size == 1)
			{
				return data;
			}

			Serializer send;
			Serializer recv;
			int mySize;
			if(rank == root)
			{
				size_t totalSize = data.size();
				size_t idealSize = totalSize / size;
				std::vector<int> counts(size, 0);
				std::vector<int> offsets(size, 0);
				size_t current = 0;
				for(rank_t _rank = 0; _rank < size; _rank++)
				{
					size_t _begin = _rank * idealSize;
					size_t _end = (_rank == size - 1)? totalSize : ((_rank + 1) * idealSize);
					size_t length = _end - _begin;
					offsets[_rank] = current;
					counts[_rank] = send.insert_elements(data, _begin, length);
					current += counts[_rank];
				}
				MPI_Scatter(counts.data(), 1, MPI_INT, &mySize, 1, MPI_INT, root, comm);
				recv.resize(mySize);
				MPI_Scatterv(send.getData(), counts.data(), offsets.data(), MPI_BYTE, recv.getData(), mySize, MPI_BYTE, root, comm);
			}
			else
			{
				MPI_Scatter(NULL, 1, MPI_INT, &mySize, 1, MPI_INT, root, comm);
				recv.resize(mySize);
				MPI_Scatterv(NULL, NULL, NULL, MPI_BYTE, recv.getData(), mySize, MPI_BYTE, root, comm);
			}

			std::vector<T> toReturn;
			recv.extract_all(toReturn);
			return toReturn;
		}

		template<typename T, typename Index_T = size_t>
		std::vector<std::vector<T>> MPI_exchange_data_indexed(const std::vector<rank_t>& correspondents, const std::vector<T>& data, const std::vector<std::vector<Index_T>> &indices = std::vector<std::vector<Index_T>>(), const size_t &extent = 1)
		{
			std::vector<MPI_Request> req(correspondents.size());
			std::vector<Serializer> senders(correspondents.size());
			for(size_t i = 0; i < correspondents.size(); ++i)
			{
				senders[i].insert_all_indexed(data, indices[i], extent);
				MPI_Isend((senders[i].size() > 0)? senders[i].getData() : NULL, senders[i].size(), MPI_CHAR, correspondents[i], MPI_EXCHANGE_TAG, MPI_COMM_WORLD, &req[i]);
			}

			std::vector<Serializer> receivers(correspondents.size());
			for(size_t i = 0; i < correspondents.size(); ++i)
			{
				MPI_Status status;
				MPI_Probe(MPI_ANY_SOURCE, MPI_EXCHANGE_TAG, MPI_COMM_WORLD, &status);
				int count;
				MPI_Get_count(&status, MPI_CHAR, &count);
				size_t location = std::distance(correspondents.begin(), std::find(correspondents.begin(), correspondents.end(), status.MPI_SOURCE));
				if(location >= correspondents.size())
				{
					MadVoro::Exception::MadVoroException eo("Bad location in mpi exchange");
					eo.addEntry("Location (Index)", location);
					eo.addEntry("Correspondents.size()", correspondents.size());
					throw eo;
				}
				receivers[location].resize(count);
				MPI_Recv(receivers[location].getData(), count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			std::vector<std::vector<T>> result(correspondents.size());
			for(size_t i = 0; i < correspondents.size(); ++i)
			{
				receivers[i].extract_all(result[i]);
			}
			if(not req.empty())
			{
				MPI_Waitall(static_cast<int>(correspondents.size()), &req[0], MPI_STATUSES_IGNORE);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			return result;
		}

	}
}

#endif // MADVORO_WITH_MPI

#endif // MPI_SERIALIZABLE_COMMANDS_HPP