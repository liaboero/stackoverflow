import sys, time, json
sys.path.insert(0, '.'); sys.path.insert(0, 'agent')
from sailing_challenge import run_episode
from forecast_agent import ForecastPlanner
params = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
for s in [7, 17, 3, 30, 40]:
    t0 = time.time()
    r = run_episode(ForecastPlanner(params), s)
    print(s, r['arrival_time'] and round(r['arrival_time'], 2), 'crashes', r['crashes'], f'{time.time()-t0:.1f}s', flush=True)
