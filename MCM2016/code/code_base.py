from __future__ import division
import numpy as np
import scipy # use numpy if scipy unavailable
import scipy.linalg # use numpy if scipy unavailable
import random
import math
import time
import numpy as np
import matplotlib
import matplotlib.animation as animation
import matplotlib.lines as mlines
matplotlib.use('TkAgg') # do this before importing pylab
import matplotlib.pyplot as plt 
import json
from pprint import pprint
import copy
#http://stackoverflow.com/questions/21352580/matplotlib-plotting-numerous-disconnected-line-segments-with-different-colors

I_factor = [24, 24, 12, 4, 1]

def lon_lat_arc_distance(lat1, long1, lat2, long2):
    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
 
    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc * 6371 #km



def gauss_2d(x0, y0, n):
    mean = [x0, y0]
    cov = [[1, 0], [0, 1]]

    x,y = np.random.multivariate_normal(mean, cov, n).T
    return zip(x,y)

class Network:
    def __init__(self):
        self.nodes = {}
        self.infected_nodes_ids = []
        self.non_infected_nodes_ids = []
        self.t = 0
        self.K = []     #how much distance affects dispersion
        self.A = 0      #captivating measure of information
        self.I = 0      #importance measure
        self.P = []
        self.B_s = 0
    
    #input:
    #node_id - name of the node
    #node_edges_ids - array of vertices id's that current node links to
    def add_node(self, node):
        node_id = node.get_id()
        input_ids = node.get_input_ids()
        output_ids = node.get_output_ids()
        if node.get_infected():
            self.infected_nodes_ids.append(node_id)
        else:
            self.non_infected_nodes_ids.append(node_id)

    def init_network(self, nodes, K, A, I, P):
        self.P = copy.deepcopy(P)
        self.K = copy.deepcopy(K)      #how much distance affects dispersion
        self.A = A      #captivating measure of information
        self.I = I      #importance measure
        self.infected_nodes_ids = []
        self.non_infected_nodes_ids = []
        self.t = 0
        self.nodes = copy.deepcopy(nodes)
        for node_id, node in nodes.iteritems():
            self.add_node(node)
        print "nodes initialized"
        print "Infected Ids", len(self.infected_nodes_ids)

    #step function for network, dt - time step 
    def step(self, dt):
        self.t += dt
        print self.t
        a_map = {} #alpha hashmap of {(i,j): edge_value}
        l_map = {} #lambda hashmpa of {i: lambda_rate(t) value}

        #for every i -> n (infected -> non_infected) edge, compute a(i,j) parameter
        for j in self.non_infected_nodes_ids:
            j_input = self.nodes[j].get_input_ids()
            l_j = 0 #lambda(i)
            for i in j_input:
                if self.nodes[i].get_infected(): #only infected neighbours affect rate
                    d = lon_lat_arc_distance(self.nodes[j].get_x(), self.nodes[j].get_y(), self.nodes[i].get_x(), self.nodes[i].get_y())  #distance between two point
                    t_i_decay = self.t - self.nodes[i].get_infection_time()
                    for i in range(len(self.K)):
                        I_i = self.I / I_factor[i]
                        bias = (math.fabs(self.nodes[j].get_bias()) + 1) / (math.fabs(self.nodes[j].get_bias() - self.B_s) + 1)
                        if bias > 0.6:
                            bias = bias**4
                            #print "BIAS", bias
                        else:
                            bias = 1
                        l_j += bias * self.P[i] * self.A*math.exp(- (t_i_decay * 1.0 / I_i + self.K[i]*d))
                               
            #print l_j
            
            if self.was_infected(l_j, dt): #infect j if exponential dist-n returned True for given dt and lambda
                self.nodes[j].infect(self.t)
                self.infected_nodes_ids.append(j)
                self.non_infected_nodes_ids.remove(j)
            
    def was_infected(self, lam, dt):
        import random
        from scipy import exp
        rnd = random.random()
        if rnd <= 1-exp(-dt*lam):
            return True
        else:
            return False    
        
    def get_infection_ratio(self):
        return len(self.infected_nodes_ids)*100.0/len(self.nodes.keys()) 

    def print_network(self):
        edges = ''
        #for key, value in self.nodes.iteritems():
        #    edges += '\t\t\t{}, input_ids: {}, output_ids: {}\n'.format(key, value.get_input_ids(), value.get_output_ids())
        print " infected_nodes_ids: \t\t\t{} \n non_infected_nodes_ids: \t\t\t{} \n node_edges: \n{} \n infection_ratio: \t\t\t{}".format(self.infected_nodes_ids, self.non_infected_nodes_ids, edges, self.get_infection_ratio())
        print "#infected - {} \n#non_infected - {}".format(len(self.infected_nodes_ids), len(self.non_infected_nodes_ids))

    def get_nodes(self):
        return self.nodes
        

class Node:
    def __init__(self, id, size=1, infection_status=0, input_ids=[], output_ids=[], x=0, y=0, time_from_infection=0):
        self.size = size
        self.id = id
        self.infected = infection_status
        self.output_ids = output_ids
        self.input_ids = input_ids
        self.x = x
        self.y = y
        self.time_from_infection = 0
        self.bias = 0

    def get_bias(self):
        return self.bias

    def get_infection_time(self):
        return self.time_from_infection

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_output_ids(self):
        return self.output_ids
    
    def get_input_ids(self):
        return self.input_ids

    def get_infected(self):
        return self.infected

    def set_bias(self, bias):
        self.bias = bias

    def set_output_ids(self, output_ids):
        self.output_ids = output_ids

    def set_input_ids(self, input_ids):
        self.input_ids = input_ids

    def infect(self, current_t):
        self.infected = 1
        self.time_from_infection = current_t
    
    def get_id(self):
        return self.id

#connect node_ids randomly
def add_random_edges(nodes, n_edges):
    node_ids = nodes.keys()
    node_in_out_ids = {} # {node_id: input, output} - to ensure input_nodes, output_nodes consistency
    for node_id in node_ids:
        node_in_out_ids.setdefault(node_id, [[],[]])
    for node_id in node_ids:
        new_node_output = [int(random.random() * len(node_ids)) for i in range(int(np.random.normal(n_edges, 1, 1)))] #normal distribution of 
        #uniqueness
        for new_id in new_node_output:           
            if (new_id not in node_in_out_ids[node_id][1]) and (new_id != node_id): #if new output_id isn't already in output_ids, add it
                #adding i-j edge on both receiving and giving side
                node_in_out_ids[new_id][0].append(node_id) #receiving side  -> j
                node_in_out_ids[node_id][1].append(new_id) #giving side     i ->

    for node_id in node_ids:
        nodes[node_id].set_input_ids(node_in_out_ids[node_id][0])
        nodes[node_id].set_output_ids(node_in_out_ids[node_id][1])
    
    return nodes

#centers = {'new york': [x,y,frac, status],...}
def make_nodes(n, centers, n_edges, infection_ratio, mode):
    nodes = {}
    inf_coordinates = []
    non_inf_coordinates = []
    inf_bias_coordinates = []
    non_inf_bias_coordinates = []
    #get coordinates
    print "getting coordinates"
    for city, center in centers.iteritems():
        #TODO: ratio of n's in different locations
        new_locations = gauss_2d(center[0], center[1], int(math.ceil(n*center[2])))
        if center[3] == 1:
            inf_coordinates += new_locations
            print len(inf_coordinates)
        elif mode == 'infect_all':
            print "INFECT ALL"
            if center[4] != 0.0:

                print "!!!!!!!!!!!!!!!!!!!!!!!!!"
                inf_bias_coordinates += new_locations
            else:
                inf_coordinates += new_locations
        else:
            #print "DONT INFECT ALL"
            if center[4] != 0.0:

                print "!!!!!!!!!!!!!!!!!!!!!!!!!"
                non_inf_bias_coordinates += new_locations
            else:
                non_inf_coordinates += new_locations

    print "creating nodes"
    for i in range(len(inf_coordinates)):
        xy = inf_coordinates[i]
        if random.random() < infection_ratio:
            infection_status = 1
        else:
            infection_status = 0
        nodes[i] = Node(i, 1, infection_status, [], [], xy[0], xy[1])

    for i in range(len(non_inf_coordinates)):
        xy = non_inf_coordinates[i]
        new_id = i + len(inf_coordinates)
        nodes[new_id] = Node(new_id, 1, 0, [], [], xy[0], xy[1])

    for i in range(len(inf_bias_coordinates)):
        xy = inf_bias_coordinates[i]
        new_id = i + len(inf_coordinates) + len(non_inf_coordinates)
        if random.random() < infection_ratio:
            infection_status = 1
        else:
            infection_status = 0
        node = Node(new_id, 1, infection_status, [], [], xy[0], xy[1])
        node.set_bias(0)
        nodes[new_id] = node

    for i in range(len(non_inf_bias_coordinates)):
        xy = non_inf_bias_coordinates[i]
        new_id = i + len(inf_coordinates) + len(non_inf_coordinates) + len(inf_bias_coordinates)
        node =  Node(new_id, 1, 0, [], [], xy[0], xy[1])
        node.set_bias(0)
        nodes[new_id] = node
    #add edges
    print "TOTAL NODES", len(nodes.items())
    print "random edgess"
    new_nodes = add_random_edges(nodes, n_edges)
    return new_nodes

def make_centers(infection_cities, bias_centers = ['Los Angeles', 'San Diego','San Francisco','San Jose'], bias_mean = 1):
    with open('cities.json') as data_file:    
        data = json.load(data_file)
    pprint(data[0]['city'])
    cities = {} # {'New York': [x, y, fraction of total people]}
    total_count = 0

    for entry in data:
        total_count += int(entry['population'])

    for entry in data:
        x = float(entry['longitude'])
        y = float(entry['latitude'])
        frac = int(entry['population']) * 1.0 / total_count
        #if frac > 0.001:
        city_name = str(entry['city'])
        if city_name in infection_cities:
            inf_status = 1
        else:
            inf_status = 0
        if city_name in bias_centers:
            bias = bias_mean
        else:
            bias = 0.0
        cities[city_name] = [x, y, frac, inf_status, bias]
    return cities

def edge_variator(centers, nodes, n_edges, n_nodes, infection_probability, infection_mode, t_max, dt,  K, A, I, P_era, n):
    n_edges = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    mode = 'asdfasdfad'
    color = 'r'
    fig = plt.figure()
    plt.axis([0,t_max,0,100])
    plt.xlabel("Time (hours)")
    plt.ylabel("Percentage of People Informed")
    for i in range(len(n_edges)):
        monte_carlo(centers, nodes, n_edges[i], n_nodes, infection_probability, mode, t_max, dt,  K, A, I, P_era, n, color, i+1)
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.show()

def monte_carlo_runner(centers, nodes, n_edges, n_nodes, infection_probability, mode, t_max, dt,  K, A, I, P_era, n, color=['r','c','b','m']): #A, I, color - arrays
    mode = 'asdfasdfad'
    fig = plt.figure()
    plt.axis([0,t_max,0,100])
    plt.xlabel("Time (hours)")
    plt.ylabel("Percentage of People Informed")
    for i in range(len(I)):
        monte_carlo(centers, nodes, n_edges, n_nodes, infection_probability, mode, t_max, dt,  K, A[i], I[i], P_era, n, color[i], i+1)
    plt.legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
    plt.show()

def monte_carlo_runner4(centers, nodes, n_edges, n_nodes, infection_probability, mode, t_max, dt,  K, A, I, P_era, n, color=['r','c','b','m']): #A, I, color - arrays
    mode = 'asdfasdfad'
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    arr = [ax1, ax2, ax3, ax4]
    for i in range(len(I)):
        nodes = make_nodes(n_nodes, centers, n_edges, infection_probability, 'infeasdfct_all')
        network = Network()
        network.init_network(nodes, K, A[i], I[i], P_era)
        
        infection_growth, changed_nodes = simulate(network, t_max, dt, '4', arr[i], i+1)
    plt.show()

def monte_carlo(centers, nodes, n_edges, n_nodes, infection_probability, mode, t_max, dt,  K, A, I, P, n, color, id):
    print 'monte carlo'
    plt.axis([0,t_max,0,100])
    #plt.title("Monte Carlo Simulation")
    plt.xlabel("Time (hours)")
    plt.ylabel("Percentage of People Informed")
    nodes = copy.deepcopy(nodes)
    for i in range(n):
        print 'monte carlo', i
        nodes = make_nodes(n_nodes, centers, n_edges, infection_probability, mode)
        network = Network()
        network.init_network(nodes, K, A, I, P)
        infection_growth, changed_nodes = simulate(network, t_max, dt, 'dont plot')
        plt.plot(np.arange(0, t_max, dt), infection_growth, linewidth=int(1/10*math.log(id)), color=color)
        if i == n-1:
            plt.plot(np.arange(0, t_max, dt), infection_growth, linewidth=int(math.log(id)), color=color, label = id)
    #plt.show()int(math.log(id))

def simulate(network, time_length, dt, mode, axis = None, id=0):
    x_all = []
    y_all = []
    color_all = []
    infection_ratio_all = []

    for i in range(int(time_length/dt)):
        changed_nodes = network.get_nodes()
        network.step(dt)

        x = []
        y = []
        color = []
        infection_ratio = network.get_infection_ratio()
        for node_id, node in changed_nodes.iteritems():
            x.append(node.get_x())
            y.append(node.get_y())
            if node.get_infected() == 1:
                color.append('r')
            else:
                color.append('b')
        x_all.append(x)
        y_all.append(y)
        color_all.append(color)
        infection_ratio_all.append(infection_ratio)
    
    if (mode == 'plot'):
        fig = plt.figure()
        #im = plt.imread('states.png')
        #implot = plt.imshow(im)
        scat = plt.scatter(x_all[0], y_all[0], c=color_all[0], s=20, alpha=0.4, label='groups')
        plt.ylabel(r'Latitude ($^{\circ}$N)', fontsize = 16)
        plt.xlabel(r'Longitude ($^{\circ}$W)', fontsize = 16)
        plt.legend(bbox_to_anchor=(0.8, 1), loc=2, borderaxespad=0.)
        colors = np.array(color_all)
        anim = animation.FuncAnimation(fig, update_plot, frames=1000,
                                      fargs=(colors, scat))
        plt.show()

        
    elif (mode == 'single'): 
        fig = plt.figure()
        #im = plt.imread('states.png')
        #implot = plt.imshow(im)
        plt.scatter(x_all[0], y_all[0], c=color_all[len(color_all)-1], s=20, alpha=0.4)
        plt.scatter(x_all[0][1], y_all[0][1], c='r', s=20, alpha=0.4, label='Informed Groups')
        plt.scatter(x_all[0][2], y_all[0][2], c='b', s=20, alpha=0.4, label='Uninformed Groups')
        print y_all[0][2]
        plt.ylabel(r'Latitude ($^{\circ}$N)', fontsize = 16)
        plt.xlabel(r'Longitude ($^{\circ}$W)', fontsize = 16)
        plt.legend(bbox_to_anchor=(0.80, 1), loc=2, borderaxespad=0.)
        colors = np.array(color_all)
        plt.show()

    elif (mode == '4'): 
        #im = plt.imread('states.png')
        #implot = plt.imshow(im)
        axis.scatter(x_all[0], y_all[0], c=color_all[len(color_all)-1], s=20, alpha=0.4)
        colors = np.array(color_all)
        axis.set_title(id, fontsize = 16)

    changed_nodes = network.get_nodes()
    return infection_ratio_all, changed_nodes

def update_plot(i, data, scat):
    scat.set_array(data[i])
    return scat,


def main():
    #dummy nodes
    dt = 2
    t_max = 96
    K = [0.67, 0.4, 0.01, 0.0001, 0.00001]
    P = [[0.75, 0, 0, 0, 0], [0.95, 0.6, 0, 0, 0], [0.85, 0.9, 0.95, 0, 0], [0.64, 0.8, 0.98, 0.14, 0], [0.3, 0.9, 0.97, 0.87, 0.64], [0.3, 0.9, 0.97, 0.99, 0.99]] # era's availability parameters, pass P[era]
    era = 2 #0 - 1870, 1 - 1920, 2- 1970, 3 - 1990, 4 - 2010
    A = 5
    I = 7
    infection_probability = 1
    n_edges = 7
    n_nodes = 1000
    n_simulations = 20
    infection_mode = 'infect_all' #infect all - makes everything infected with probability infection_probability, else - only infects infection_cities
    infection_cities = ['Los Angeles']
    print "importing locations"
    centers = make_centers(infection_cities)
    nodes = make_nodes(n_nodes, centers, n_edges, infection_probability, infection_mode) #bias in here
    network = Network()
    network.init_network(nodes, K, A, I, P[era])
    
    
    infection_growth, changed_nodes = simulate(network, t_max, dt, 'single') #bias in here
    plt.axis([0,t_max,0,100])
    plt.plot(np.arange(0, t_max, dt), infection_growth, linewidth=1, color='r')
    plt.show()
    print "simulating spread"
    n_nodes = 1000


    A_array0 = [10, 10, 40, 1000]       #dt = 2, t_max = 96, i = 1
    I_array0 = [15, 50, 15, 500]
    A_array1 = [100, 100, 500, 1000]    #dt = 2, t_max = 96, i = 1
    I_array1 = [100, 500, 100, 500]
    A_array2 = [10, 10, 40, 40]         #dt = 2, t_max = 96, i = 1
    I_array2 = [15, 50, 15, 50]
    A_array3 = [2, 2, 5, 5]             #ip = 0.1
    I_array3 = [2, 7, 2, 7]
    A_array4 = [0.1, 0.1, 1, 1]         #ip = 1
    I_array4 = [0.3, 1, 0.3, 1]
    #monte_carlo_runner(centers, nodes, n_edges, n_nodes, infection_probability, infection_mode, t_max, dt,  K, A_array4, I_array4, P[era+1], n_simulations)
    n_nodes = 5000   
    monte_carlo_runner4(centers, nodes, n_edges, n_nodes, infection_probability, infection_mode, t_max, dt,  K, A_array2, I_array2, P[era], n_simulations)
    #edge_variator(centers, nodes, n_edges, n_nodes, infection_probability, infection_mode, t_max, dt,  K, 1, 1, P[era], 1)
    


main()