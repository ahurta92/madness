/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

//#define MAD_ARCHIVE_DEBUG_ENABLE

#include <algorithm>

#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <madness/world/worldmutex.h>
#include <madness/world/atomicint.h>

#include <madness/world/vector_archive.h>
#include <madness/world/parallel_archive.h>

using namespace madness;
using namespace std;

struct Key {
    int k;

    Key() : k(-1) {}

    Key(int k) : k(k) {}

    hashT hash() const {
        return k;
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & k;
    }

    bool operator==(const Key& b) const {
        return k==b.k;
    }
};

ostream& operator<<(ostream&s, const Key& key) {
    s << "Key(" << key.k << ")";
    return s;
}

struct Node {
    int k;

    Node() : k(-1) {}

    Node(int k) : k(k) {}

    int get() const {
        return k;
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & k;
    }

    ~Node() {}
};

ostream& operator<<(ostream&s, const Node& node) {
    s << "Node(" << node.k << ")";
    return s;
}


void test0(World& world) {
    WorldContainer<Key,Node> c(world);

    Key key1(1);
    Node node1(1);

    if (c.owner(key1) == world.rank()) c.replace(key1,node1);

    world.gop.fence();

    for (int i=0; i<10000; ++i)
        MADNESS_CHECK(c.find(key1).get()->second.get() == 1);

    for (int i=3; i<100; ++i)
        MADNESS_CHECK(c.find(Key(i)).get() == c.end());

    world.gop.fence();
}

class TestPmap : public WorldDCPmapInterface<int> {
private:
    const int nproc;
    const int shift;
public:
    TestPmap(World& world, int shift)
        : nproc(world.mpi.nproc())
        , shift(shift)
    { }

    ProcessID owner(const int& key) const {
        if (nproc == 1) return 0;
        return (key + shift)%nproc;
    }
};

AtomicInt double_count;

// Can use just a plain int as key or value but make this class only to conunt instances
class Double {
    double value;
public:

    Double(double value=0.0) : value(value) {double_count++;}

    Double(const Double& d) : value(d.value) {double_count++;}

    ~Double() {double_count--;}

    Double operator+(double x) {return Double(value+x);}

    bool operator==(double x) {return value==x;}


    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & value;
    }
};
    

void test1(World& world) {
    std::shared_ptr< WorldDCPmapInterface<int> > pmap0(new TestPmap(world, 0));
    std::shared_ptr< WorldDCPmapInterface<int> > pmap1(new TestPmap(world, 1));

    int total;
    double_count = 0;
    world.gop.fence();
        
    {
        WorldContainer<int,Double> c(world,pmap0), d(world,pmap0), e(world,pmap0);
        
        world.gop.fence(); total = double_count; world.gop.sum(total);
        if (world.rank() == 0) print("count after constructor", total);
        
        if (world.rank() == 0) {
            for (int i=0; i<100; ++i) {
                c.replace(i,i+1.0);
                d.replace(i,i+2.0);
                e.replace(i,i+3.0);
            }
        }
        
        world.gop.fence(); total = double_count; world.gop.sum(total);
        if (world.rank() == 0) print("count after making", total);
        
        pmap0->redistribute(world, pmap1);
        
        world.gop.fence(); total = double_count; world.gop.sum(total);
        if (world.rank() == 0) print("count after redistributing", total);
        std::size_t global_count = c.get_pmap()->global_size(world);
        if (world.rank() == 0) print("count after from global sz", global_count);
        
        for (int i=0; i<100; ++i) {
            MADNESS_CHECK(c.find(i).get()->second == (i+1.0));
            MADNESS_CHECK(d.find(i).get()->second == (i+2.0));
            MADNESS_CHECK(e.find(i).get()->second == (i+3.0));
        }
        
        world.gop.fence(); total = double_count; world.gop.sum(total);
        if (world.rank() == 0) print("count after testing", total);
        
        c.clear();
        
        world.gop.fence(); total = double_count; world.gop.sum(total);
        if (world.rank() == 0) print("count after clearing c", total);
    }

    total = double_count; world.gop.sum(total);
    if (world.rank() == 0) print("count after going out of scope without fence", total);
    
    world.gop.fence(); total = double_count; world.gop.sum(total);
    if (world.rank() == 0) print("count after first fence", total);

    world.gop.fence(); total = double_count; world.gop.sum(total);
    if (world.rank() == 0) print("count after second fence", total);
}


void test_local(World& world) {

	print("entering test_local");
    std::shared_ptr< WorldDCPmapInterface<int> > pmap0(new TestPmap(world, 0));

    WorldContainer<int,double> c(world,pmap0),d(world,pmap0),e(world,pmap0);
    for (int i=0; i<10; ++i) {
    	c.replace(i,i+1.0);
    	d.replace(i,i+1.0);
    	e.replace(i,i+1.0);
    }
    world.gop.fence();

    std::size_t size=c.size();
    if (world.rank()==0) print("size on rank=0 before replication",size);

    c.replicate();		// fence
    d.replicate(false);	// no fence
    e.replicate();		// fence

    size=c.size();
    if (world.rank()==0) print("size on rank=0 after replication",size);

    auto localmap=c.get_pmap();
    localmap->redistribute(world, pmap0);
    world.gop.fence();
    size=c.size();
    if (world.rank()==0) print("size on rank=0 after redistribution",size);

    localmap=d.get_pmap();
    localmap->redistribute(world, pmap0);
    localmap->redistribute(world, pmap0);

}

namespace madness {
    namespace archive {
        /// Write container to parallel archive
        template <class keyT, class valueT>
        struct ArchiveStoreImpl< ParallelOutputArchive<VectorOutputArchive>, WorldContainer<keyT,valueT> > {
            static void store(const ParallelOutputArchive<VectorOutputArchive>& ar, const WorldContainer<keyT,valueT>& t) {
                using localarchiveT = VectorOutputArchive;
                const long magic = -5881828; // Sitar Indian restaurant in Knoxville (negative to indicate parallel!)
                typedef WorldContainer<keyT,valueT> dcT;
                using const_iterator = typename dcT::const_iterator;

                const size_t default_size = 100*1024*1024;

                World* world = ar.get_world();
                world->gop.fence();

                std::vector<unsigned char> v;
                v.reserve(default_size);
                size_t count = 0;

                class op : public TaskInterface {
                    const size_t ntasks;
                    const size_t taskid;
                    const dcT& t;
                    std::vector<unsigned char>& vtotal;
                    size_t& total_count;
                    Mutex& mutex;
                    
                public:
                    op(size_t ntasks, size_t taskid, const dcT& t, std::vector<unsigned char>& vtotal, size_t& total_count, Mutex& mutex)
                        : ntasks(ntasks), taskid(taskid), t(t), vtotal(vtotal), total_count(total_count), mutex(mutex) {}
                    void run(World& world) {
                        std::vector<unsigned char> v;
                        v.reserve(std::max(size_t(1024*1024),vtotal.capacity()/ntasks));
                        VectorOutputArchive var(v);
                        const_iterator it=t.begin();
                        size_t count = 0;
                        size_t n = 0;
                        while (it!=t.end()) {
                            if ((n%ntasks) == taskid) {
                                var & *it;
                                ++count;
                            }
                            ++it;
                            n++;
                        }

                        if (count) {
                            mutex.lock();
                            vtotal.insert(vtotal.end(), v.begin(), v.end());
                            total_count += count;
                            mutex.unlock();
                        }
                    }
                };
                
                Mutex mutex;
                size_t ntasks = std::max(size_t(1), ThreadPool::size());
                for (size_t taskid=0; taskid<ntasks; taskid++)
                    world->taskq.add(new op(ntasks, taskid, t, v, count, mutex));
                world->gop.fence();

                // Gather all buffers to process 0
                // first gather all of the sizes and counts to a vector in process 0
                int size = v.size();
                std::vector<int> sizes(world->size());
                MPI_Gather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, world->mpi.comm().Get_mpi_comm());
                world->gop.sum(count); // just need total number of elements

                // build the cumulative sum of sizes
                std::vector<int> offsets(world->size());
                offsets[0] = 0;
                for (int i=1; i<world->size(); ++i) offsets[i] = offsets[i-1] + sizes[i-1];
                int total_size = offsets.back() + sizes.back();

                // gather the vector of data v from each process to process 0
                unsigned char* all_data=0;
                if (world->rank() == 0) {
                    all_data = new unsigned char[total_size];
                }
                MPI_Gatherv(v.data(), v.size(), MPI_BYTE, all_data, sizes.data(), offsets.data(), MPI_BYTE, 0, world->mpi.comm().Get_mpi_comm());

                if (world->rank() == 0) {
                    auto& localar = ar.local_archive();
                    localar & magic & 1; // 1 client
                    // localar & t;
                    ArchivePrePostImpl<localarchiveT,dcT>::preamble_store(localar);
                    localar & -magic & count;                
                    localar.store(all_data, total_size);
                    ArchivePrePostImpl<localarchiveT,dcT>::postamble_store(localar);

                    delete[] all_data;
                }
                world->gop.fence();
            }
        };
    }
}

void test_florian(World& world) {
    WorldContainer<Key,Node> c(world);

    Key key1(1);
    Node node1(1);

    if (world.rank() == 0) {
        for (int i=0; i<100; ++i) {
            c.replace(Key(i),Node(i));
        }
    }
    world.gop.fence();

    std::vector<unsigned char> v;
    {
        archive::VectorOutputArchive var(v);
        archive::ParallelOutputArchive ar(world,var);
        ar & c;
    }

    WorldContainer<Key,Node> c2(world);
    {
        archive::VectorInputArchive var2(v);
        archive::ParallelInputArchive ar2(world,var2);
        ar2 & c2;
    }

    for (int i=0; i<100; ++i) {
        MADNESS_CHECK(c2.find(Key(i)).get()->second.get() == i);
    }

    world.gop.fence();
    print("test_florian passed");
}

int main(int argc, char** argv) {

    try {
        World& world = initialize(argc, argv);
        // test0(world);
        // test1(world);
        // test1(world);
        // test1(world);
        // test_local(world);
        test_florian(world);
    }
    catch (const SafeMPI::Exception& e) {
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a string exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    finalize();
    return 0;
}
