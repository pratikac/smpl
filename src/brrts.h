#ifndef __brrts_h__
#define __brrts_h__

#include <list>
#include <set>
#include <vector>
#include <set>
#include <cfloat>
#include "kdtree.h"
#include <algorithm>
#include <tuple>
#include <functional>
#include <unordered_map>

#include "system.h"
#include "utils.h"

#include <lcm/lcm.h>
#include <bot_core/bot_core.h>
#include <bot_lcmgl_client/lcmgl.h>
#include <bot_vis/gl_util.h>

using namespace std;

template<class system_tt> class bvertex_c;
template<class system_tt> class bedge_c;

template<class system_tt>
class bvertex_c
{
  public:
    typedef bvertex_c<system_tt> bvertex;
    typedef bedge_c<system_tt> bedge;

    typedef system_tt  system_t;
    typedef typename system_tt::state state_t;
    typedef typename system_tt::control control_t;
    typedef typename system_tt::opt_data_t opt_data_t;
    typedef typename system_tt::cost_t cost_t;
    typedef typename system_tt::trajectory trajectory_t;

    bvertex* child;
    state_t state;
    double t0;
    set<bvertex*> parents;

    int mark;

    cost_t cost_to_root;
    cost_t cost_to_child;
    bedge* bedge_to_child;

    bvertex_c()
    {
      child = NULL;
      bedge_to_child = NULL;
      mark = 0;
      t0 =0;
    }
    ~bvertex_c()
    {
      if(bedge_to_child)
        delete bedge_to_child;
    }
    bvertex_c(const state_t& si)
    {
      child = NULL;
      bedge_to_child = NULL;
      state = si;
      mark = 0;
      t0 =0;
    }
    
    string tostring() const
    {
      return state->tostring();
    }
    void print_branch(string prefix)
    {
      cout << prefix << tostring() << " - cost: "; cost_to_root->print();
      for(auto& pc : parents)
        pc->print_branch(prefix+"\t");
    }
    state_t& get_state() const { return state;};
    bvertex& get_child() const {return *child;};
    cost_t& get_cost() const {return cost_to_root;};
};

template<class system_tt>
class bedge_c
{
  public:
    typedef bedge_c<system_tt> bedge;

    typedef system_tt  system_t;
    typedef typename system_tt::state state;
    typedef typename system_tt::control control;
    typedef typename system_tt::opt_data_t opt_data_t;
    typedef typename system_tt::trajectory trajectory;
    typedef typename system_tt::cost_t cost_t;  

    const state* start_state;
    const state* end_state;

    cost_t cost;
    opt_data_t opt_data;
    double dt;

    bedge_c()
    {
      end_state = NULL;
      start_state = NULL;
      dt = 0;
    };

    bedge_c(const state* si, const state* se, cost_t& c, double dt_in, opt_data_t& opt_data_in)
    {
      start_state = si;
      end_state = se;
      opt_data = opt_data_in;
      cost = c;
      dt = dt_in;
    }
};

template<class bvertex_tt, class bedge_tt>
class brrts_c
{
  public:
    typedef typename bvertex_tt::system_t system_t;
    typedef typename system_t::state state;
    typedef typename system_t::control control;
    typedef typename system_t::opt_data_t opt_data_t;
    typedef typename system_t::trajectory trajectory_t;
    
    typedef typename system_t::cost_t cost_t;  
    typedef typename system_t::region_t region_t;

    const static size_t num_dim = system_t::N;

    typedef bvertex_tt bvertex;
    typedef bedge_tt bedge;

    typedef struct kdtree kdtree_t;
    typedef struct kdres kdres_t;

    system_t system;

    int num_vertices;
    list<bvertex*> list_vertices;
    double gamma;
    double goal_sample_freq;
    bool do_branch_and_bound;

    bvertex* root;
    cost_t lower_bound_cost;
    bvertex* lower_bound_bvertex;
    kdtree_t* kdtree;
    bvertex* last_added_bvertex;

    static int debug_counter;
    bot_lcmgl_t* lcmgl;
    double points_color[4];
    double points_size;
    double lines_color[4];
    double lines_width;
    double best_lines_color[4];
    double best_lines_width;

    brrts_c(){
      basic_initialization();
    }
   
    void set_points_color(double* pc, double ps)
    {
      for(int i : range(0,4))
        points_color[i] = pc[i];
      points_size = ps;
    }
    
    void set_lines_color(double* lc, double lw)
    {
      for(int i : range(0,4))
        lines_color[i] = lc[i];
      lines_width = lw;
    }

    void set_best_lines_color(double* blc, double blw)
    {
      for(int i : range(0,4))
        best_lines_color[i] = blc[i];
      best_lines_width = blw;
    }
    
    void basic_initialization()
    {
      gamma = 2.5;
      goal_sample_freq = 0.1;
      do_branch_and_bound = true;

      root = NULL;
      lower_bound_bvertex = NULL;
      last_added_bvertex = NULL;

      kdtree = NULL;
      num_vertices = 0;
  
      double pc[4] = {1,1,0,0.9};
      double lc[4] = {1,1,1,0.5};
      double blc[4] = {1,0,0,0.8};
      set_points_color(pc, 4);
      set_lines_color(lc, 1.5);
      set_best_lines_color(blc, 4);
    }
    
    brrts_c(bot_lcmgl_t* lcmgl_in){
      lcmgl = lcmgl_in;
      basic_initialization();
    }
    
    ~brrts_c()
    {
      if(kdtree)
        kd_free(kdtree);
      clear_list_vertices();
    }

    void clear_list_vertices()
    {
      for(auto& i : list_vertices)
        delete i;
      list_vertices.clear();
      num_vertices = 0;
    }
    
    int set_root(const state& rs)
    {
      root = new bvertex(rs);
      root->cost_to_root = system.get_zero_cost();
      root->cost_to_child = system.get_zero_cost();
      root->bedge_to_child = NULL;
      root->child = NULL;
      root->parents.clear();

      insert_into_kdtree(*root);
      //update_best_bvertex(*root);
      return 0; 
    }
    
    int initialize(const state& rs, bool do_branch_and_bound_in=true)
    {
      clear_list_vertices();  
      lower_bound_cost = system.get_inf_cost();
      lower_bound_bvertex = NULL;
      do_branch_and_bound = do_branch_and_bound_in;

      if(kdtree)
        kd_free(kdtree);
      kdtree = kd_create(num_dim);

      set_root(rs);
      root->state.print(cout,"set root to:", "\n");
      last_added_bvertex = root;

      return 0;
    }

    int check_collision_trajectory(const trajectory_t& t1, const trajectory_t& t2, double dmax)
    {
      double dmin = FLT_MAX/2;
      double t = min(t1.t0, t2.t0);
      double dt = t1.dt;
      double T = min(t1.t0 + t1.dt*t1.states.size(), t2.t0 + t2.dt*t2.states.size());
      
      int i = 0;
      while(t < T)
      {
        auto& s1 = t1.states[i];
        auto& s2 = t2.states[i];
        double t3 = s1.dist(s2);
        if(t3 < dmax)
          return true;
        
        i += 10;
        t += 10*dt;
      }
      return false;
    }
        
    int iteration(state* s_in = NULL, set<bvertex*>* rewired_vertices=NULL, trajectory_t* obstacle_trajectory=NULL, double collision_distance = 1)
    {
      last_added_bvertex = NULL;

      // 1. sample
      state sr;
      if(!s_in)
      {
        int ret = 0;
        double p = RANDF;
        if(p < goal_sample_freq)
          ret = system.sample_in_goal(sr);
        else
          ret = system.sample_state(sr);
        if(ret)
          return 1;
      }
      else
        sr = *s_in;

      // 2. compute nearest vertices
      vector<bvertex*> near_vertices;
      if(get_near_vertices(sr, near_vertices))
        return 2;

      // 3. best child
      bvertex* best_child = NULL;
      bedge* bedge_to_child = NULL;
      if(find_best_child(sr, near_vertices, best_child, bedge_to_child))
        return 3;

      // 4.a check if the trajectory new sample collides with collision_trajectory
      if(obstacle_trajectory)
      {
        trajectory_t tmp_traj;
          system.extend_to(sr, best_child->state, false,
              tmp_traj, bedge_to_child->opt_data);
        tmp_traj.t0 = best_child->t0; 

        if(check_collision_trajectory(*obstacle_trajectory, tmp_traj, collision_distance))
          return 4;
      }
      
      // 4. draw bedge to child from new bvertex
      bvertex* new_bvertex = insert_bedge(*best_child, *bedge_to_child);
      if(!new_bvertex)
        return 5;

      // 5. rewire
      if(near_vertices.size())
        rewire_vertices(*new_bvertex, near_vertices, rewired_vertices);
      
      last_added_bvertex = new_bvertex;
      return 0;
    }

    int insert_into_kdtree(bvertex& v)
    {
      double* key = new double[num_dim];
      system.get_key(v.state, key);
      kd_insert(kdtree, key, &v);
      delete[] key;

      list_vertices.push_back(&v);
      num_vertices++;
      return 0;
    }

    bvertex& get_root_bvertex() {return *root;};
    cost_t get_best_cost()    {return lower_bound_cost;};
    bvertex& get_best_bvertex() {return *lower_bound_bvertex;}

    int get_trajectory_root(bvertex& v, trajectory_t& root_traj)
    {
      root_traj.clear();
      bool check_obstacles = false;
      bvertex* vc = static_cast<bvertex*>(&v);
      while(vc)
      {
        bvertex* vchild = static_cast<bvertex*>(vc->child);
        if(vchild)
        {
          trajectory_t traj_to_child;
          system.extend_to(vc->state, vchild->state, check_obstacles,
              traj_to_child, vc->bedge_to_child->opt_data);
          root_traj.append(traj_to_child);
        }
        vc = vchild;
      }
      root_traj.t0 = 0;
      return 0;
    }

    int get_best_trajectory(trajectory_t& best_traj)
    {
      if(!lower_bound_bvertex)
        return 1;
        
      get_trajectory_root(*lower_bound_bvertex, best_traj);
      return 0;
    }

    int get_nearest_vertex(const state& s, bvertex*& nearest_vertex)
    {
      int toret = 0;
      double* key = new double[num_dim];
      system.get_key(s, key);

      double rn = gamma*pow(log(num_vertices + 1.0)/(num_vertices+1.0), 1.0/(double)num_dim);
      kdres_t* kdres = kd_nearest(kdtree, key);
      if(!kd_res_size(kdres))
        toret = 1;
      else
        nearest_vertex = (bvertex*) kd_res_item_data(kdres);
      
      delete[] key;
      kd_res_free(kdres);
      return toret;
    }
    
    int get_near_vertices(const state& s, vector<bvertex*>& near_vertices)
    {
      int toret = 0;
      double* key = new double[num_dim];
      system.get_key(s, key);

      double rn = gamma*pow(log(num_vertices + 1.0)/(num_vertices+1.0), 1.0/(double)num_dim);
      kdres_t* kdres = kd_nearest_range(kdtree, key, rn);

      int num_near_vertices = kd_res_size(kdres);
      if(!num_near_vertices)
      {
        kd_res_free(kdres);

        // get nearest bvertex
        kdres = kd_nearest(kdtree, key);
        if(kd_res_end(kdres))
          toret = 1;
        else
        {
          bvertex* vc = (bvertex*) kd_res_item_data(kdres);
          near_vertices.push_back(vc);
        }
      }
      else
      {
        while(! kd_res_end(kdres))
        {
          bvertex* vc = (bvertex*)kd_res_item_data(kdres);
          near_vertices.push_back(vc);
          kd_res_next(kdres);
        }
      }

      delete[] key;
      kd_res_free(kdres);
      return toret;
    }

    int update_best_bvertex(bvertex& v)
    {
      if(system.is_in_goal(v.state))
      {
        // implement goal cost here
        if( (!lower_bound_bvertex) || (v.cost_to_root < lower_bound_cost))
        {
          lower_bound_cost = v.cost_to_root;
          lower_bound_bvertex = &v;
        }
      }
      return 0;
    }

    bvertex* insert_bedge(bvertex& ve, bedge& e)
    {
      // branch and bound
      if(do_branch_and_bound)
      {
        cost_t new_cost = ve.cost_to_root + e.cost; 
        if(new_cost > lower_bound_cost)
          return NULL;
      }

      // create new bvertex
      bvertex* new_bvertex = new bvertex(*(e.start_state));
      insert_into_kdtree(*new_bvertex);

      insert_bedge(*new_bvertex, e, ve);
      return new_bvertex;
    }

    int insert_bedge(bvertex& vp, bedge& e, bvertex& vc)
    {
      vc.t0 = vp.t0 + e.dt;

      vp.cost_to_child = e.cost;
      vp.cost_to_root = vc.cost_to_root + vp.cost_to_child;
      update_best_bvertex(vp);

      if(vp.bedge_to_child)
        delete vp.bedge_to_child;
      vp.bedge_to_child = &e;

      if(vp.child)
        vp.child->parents.erase(&vp);
      vp.child = &vc;
      vc.parents.insert(&vp);
      return 0;
    }

    static bool compare_bvertex_cost_pairs(const pair<bvertex*, cost_t>& p1,
        const pair<bvertex*, cost_t>& p2)
    {
      return (p1.second < p2.second);
    }

    int find_best_child(const state& si, const vector<bvertex*>& near_vertices,
        bvertex*& best_child, bedge*& best_bedge)
    {
      // 1. create bvertex_cost_pairs
      vector<pair<bvertex*, cost_t> > bvertex_cost_pairs;
      unordered_map<bvertex*, tuple<cost_t, cost_t, opt_data_t> > bvertex_map;
      for(auto& pv : near_vertices)
      {
        bvertex& v = *pv;
        opt_data_t opt_data;
        cost_t bedge_cost;

        if(system.evaluate_extend_cost(si, v.state, opt_data, bedge_cost))
          continue;
        cost_t v_cost = v.cost_to_root + bedge_cost;

        bvertex_cost_pairs.push_back(make_pair(pv, v_cost));
        bvertex_map.insert(make_pair(pv, make_tuple(bedge_cost, v_cost, opt_data)));
      }

      // 2. sort using compare function of cost_t
      sort(bvertex_cost_pairs.begin(), bvertex_cost_pairs.end(), compare_bvertex_cost_pairs);

      // 3. check obstacles in order of increasing cost
      bool check_obstacles = true;
      trajectory_t traj;
      for(auto& p : bvertex_cost_pairs)
      {
        bvertex& v = *(p.first);
        opt_data_t& opt_data = get<2>(bvertex_map[p.first]);
        cost_t& bedge_cost = get<0>(bvertex_map[p.first]);
        if(!system.extend_to(si, v.state, check_obstacles,traj, opt_data))
        {
          double edge_duration = system.dynamical_system.evaluate_extend_cost(si, v.state, opt_data);
          best_child = &v;
          best_bedge = new bedge(&si, &(v.state), bedge_cost, edge_duration, opt_data);
          //cout<<"best_bedge.cost: "<< best_bedge.cost.val << endl;
          return 0;
        }
      }
      return 1;
    }

    int update_all_costs()
    {
      lower_bound_cost = system.get_inf_cost();
      lower_bound_bvertex = NULL;
      update_branch_cost(*root,0); 
      return 0;
    }

    int update_branch_cost(const bvertex& v, int depth)
    {
      for(auto& pc : v.parents)
      {
        bvertex& parent = *(static_cast<bvertex*>(pc));
        parent.cost_to_root = v.cost_to_root + parent.cost_to_child;
        update_best_bvertex(parent);
        update_branch_cost(parent, depth+1);
      }
      return 0;
    }

    int rewire_vertices(bvertex& v, const vector<bvertex*>& near_vertices, set<bvertex*>* rewired_vertices)
    {
      bool check_obstacles = true;
      for(auto& pvn : near_vertices)
      {
        bvertex& vn = *pvn;
        opt_data_t opt_data;
        trajectory_t traj;
        cost_t cost_bedge;
        if(system.evaluate_extend_cost(vn.state, v.state, opt_data, cost_bedge))
          continue;

        if(rewired_vertices)
          rewired_vertices->insert(pvn);
        cost_t cvn = v.cost_to_root + cost_bedge;
        if(cvn < vn.cost_to_root)
        {
          if(system.extend_to(vn.state, v.state, check_obstacles, traj, opt_data))
            continue;
          
          double en_dt = system.dynamical_system.evaluate_extend_cost(vn.state, v.state, opt_data);
          bedge* en = new bedge(&(vn.state), &(v.state), cost_bedge, en_dt, opt_data);
          insert_bedge(vn, *en, v);

          update_branch_cost(vn, 0);
        }
      }
      return 0;
    }

    int recompute_cost(bvertex& v)
    {
      update_branch_cost(v,0);  
      return 0;
    }

    bool is_safe_trajectory(const trajectory_t& traj)
    {
      return system.is_safe_trajectory(traj);
    }

    int mark_ancestor_vertices(bvertex& v)
    {
      v.mark = 1;
      for(auto& pc : v.parents)
        mark_ancestor_vertices(*pc);
      return 0;
    }
    
    int delete_unmarked_vertices(vector<bvertex*>& surviving_vertices)
    {
      surviving_vertices.clear();
      for(auto& pv : list_vertices)
      {
        if(pv->mark == 1)
        {
          pv->mark = 0;
          surviving_vertices.push_back(pv);
        }
        else
          delete pv;
      }
      return 0;
    }
    
    int mark_bvertex_and_remove_from_child(bvertex& v)
    {
      v.mark = 1;
      v.child->parents.erase(&v);
      for(auto& pc : v.parents)
        mark_bvertex_and_remove_from_child(*pc);
      return 0;
    }

    int check_and_mark_parents(bvertex& v)
    {
      trajectory_t traj;
      bool check_obstacles = true;

      if(system.extend_to(v.child->state, v.state, check_obstacles,
            traj, v.bedge_to_child->opt_data))
      { 
        mark_bvertex_and_remove_from_child(v);
      }
      else
      {
        for(auto& pc : v.parents)
          check_and_mark_parents(*pc);
      }
      return 0;
    }

    int check_tree()
    {
      if(system.is_in_collision(root->state))
      {
        cout<<"root in collision"<<endl;
        return 1;
      }
      else if(!root->parents.empty())
      {
        root->mark = 0;

        for(auto& prc : root->parents)
          check_and_mark_parents(*prc);

        list<bvertex*> surviving_vertices;
        for(auto& pv : list_vertices)
        {
          if(!pv->mark)
            surviving_vertices.push_back(pv);
          else
            delete pv;
        }

        if(kdtree)
          kd_free(kdtree);
        kdtree = kd_create(num_dim);

        list_vertices.clear();
        num_vertices = 0;
        for(auto& pv : surviving_vertices)
          insert_into_kdtree(*pv); 

        update_all_costs();
      }
      return 0;
    }

    int lazy_check_tree(const trajectory_t& committed_trajectory)
    {
      int ret = 0;
      if(!system.is_safe_trajectory(committed_trajectory))
        ret = check_tree();
      return ret;
    }

    int get_best_trajectory_vertices(list<bvertex*>& best_trajectory_vertices)
    {
      if(!lower_bound_bvertex)
        return 1;
      best_trajectory_vertices.clear();
      bvertex* pvc = lower_bound_bvertex;
      while(pvc)
      {
        best_trajectory_vertices.push_back(pvc);
        pvc = pvc->child;
      }
      return 0;
    }

    int print_marks()
    {
      cout<<"marks: ";
      for(auto& pv : list_vertices)
        cout<<pv->mark<<" ";
      cout<<endl;
      return 0;
    }
    
    virtual void plot_trajectory(trajectory_t& traj, double* lc, double width)
    {
      bot_lcmgl_color4f(lcmgl, lc[0], lc[1], lc[2], lc[3]);
      bot_lcmgl_line_width(lcmgl, width);
      bot_lcmgl_begin(lcmgl, GL_LINES);
      for(auto pts = traj.states.begin(); pts != traj.states.end(); pts++)
      {
        double s2[3] = {0};
        system.get_plotter_state(*pts, s2);
        bot_lcmgl_vertex3d(lcmgl, s2[0], s2[1], s2[2]);

        auto npts = pts;
        npts++;
        if(npts != traj.states.end()){
          system.get_plotter_state(*npts, s2);
          bot_lcmgl_vertex3d(lcmgl, s2[0], s2[1], s2[2]);
        }
      }
      bot_lcmgl_end(lcmgl);
    }

    virtual void plot_tree()
    {
      bool check_obstacles = false;
      if(num_vertices == 0)
        return;
      for(auto& v : list_vertices)
      {
        double s1[3] = {0};
        system.get_plotter_state(v->state, s1);
        bot_lcmgl_color4f(lcmgl, points_color[0], points_color[1], points_color[2], points_color[3]);
        bot_lcmgl_point_size(lcmgl, points_size);
        bot_lcmgl_begin(lcmgl, GL_POINTS);
        bot_lcmgl_vertex3d(lcmgl, s1[0], s1[1], s1[2]);
        bot_lcmgl_end(lcmgl);
        
        if(v->child){
          trajectory_t traj_to_child;
          if(system.extend_to(v->state, v->child->state, check_obstacles,
              traj_to_child, v->bedge_to_child->opt_data)){
            cout<<"extend_to returns 1 while plotting"<<endl;
            return;
          }
          plot_trajectory(traj_to_child, lines_color, lines_width);
        }
      }
    }

    virtual void plot_best_trajectory()
    {
      trajectory_t best_trajectory;
      if(get_best_trajectory(best_trajectory))
        return;
      plot_trajectory(best_trajectory, best_lines_color, best_lines_width); 
    }

    void plot_region(region_t& r)
    {
      double c[3] = {0};
      double s[3] = {0};
      float sf[3] = {0};
      r.get_plotter_state(c, s);
      bot_lcmgl_color4f(lcmgl, r.color[0], r.color[1], r.color[2], r.color[3]);
      double cd[3] = {0};
      for(int i=0; i<3; i++)
        sf[i] = s[i];
      bot_lcmgl_box(lcmgl, c, sf);
    }

    virtual void plot_environment()
    {
      bot_lcmgl_enable(lcmgl, GL_BLEND);

      plot_region(system.operating_region);
      plot_region(system.goal_region);
    }
};
#endif
