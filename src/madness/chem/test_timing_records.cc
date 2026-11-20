#include <madness/chem/TimingRecords.h>
#include <stdexcept>
using namespace madness::iteration;
struct ManualClock {
    using duration = std::chrono::milliseconds;
    using time_point = std::chrono::time_point<ManualClock>;
    static inline time_point current{};
    static time_point now() { return current; }
    static void advance(int milliseconds) { current += duration(milliseconds); }
};
void require(bool ok) { if (!ok) throw std::runtime_error("timing record contract failed"); }
int main() {
    TimingCollector<ManualClock> timing;
    timing.begin("run", {{"clock_scope", "driver_rank"}});
    ManualClock::advance(100);
    timing.begin("response", {{"calculation_id", "response"}});
    timing.begin("task", {{"task_id", "coupled"}, {"execution", {{"mpi_ranks", 2}, {"threads_per_rank", 3}}}});
    timing.begin("evaluate", {{"iteration", 0}});
    ManualClock::advance(200);
    timing.end();
    timing.unavailable("bsh", "included_in_evaluate");
    ManualClock::advance(300);
    timing.end();
    timing.end();
    ManualClock::advance(400);
    const auto live = timing.record();
    require(live.at("seconds") == 1.0);
    timing.end();
    require(timing.record() == live);
    const auto& task = live.at("children").at(0).at("children").at(0);
    require(task.at("seconds") == 0.5);
    require(task.at("execution").at("mpi_ranks") == 2);
    require(task.at("children").at(0).at("seconds") == 0.2);
    require(!task.at("children").at(1).contains("seconds"));
    require(task.at("children").at(1).at("availability").at("status") == "unavailable");
    TimingCollector<ManualClock> restarted;
    restarted.begin("run", {{"clock_scope", "driver_rank"}});
    restarted.begin("load", {{"task_id", "coupled"}});
    ManualClock::advance(25);
    restarted.end(); restarted.end();
    require(restarted.record().at("seconds") == 0.025);
    require(restarted.record().at("children").size() == 1);
    bool rejected = false;
    try { restarted.end(); } catch (const std::exception&) { rejected = true; }
    require(rejected);
}
