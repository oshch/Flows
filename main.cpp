#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>
#include <stack>
#include <list>

using namespace std;

const int INF = 1e9;
const long long INFLL = (long long)1e16;

class Edge {
protected:
	size_t from, to;
	size_t id;
	bool used;
public:
	Edge() {}
	Edge(size_t _from, size_t _to, size_t _id) : from(_from), to(_to), id(_id) {}
	void setUsed(bool u);
	bool getUsed();
	size_t getTo();
	size_t getId();
	size_t getFrom();
};

class Vertex {
protected:
	size_t dist, id;
	list<Edge*>::iterator edgesHead;
	long long excess;
	int label;
public:
	Vertex() : dist(INF), edgesHead(NULL) {}

	size_t getDist();
	list<Edge*>::iterator getHead();
	void setDist(size_t x);
	void increaseHead();
	void initHead(list<Edge*>::iterator i);

	long long getExcess();
	int getLabel();
	void setExcess(long long x);
	void setLabel(int x);
	void increaseExcess(long long x);
	void setId(size_t x);
	size_t getId();
};

class FlowEdge : public Edge {
protected:
	int flow, capacity;
	FlowEdge* rev;
public:
	FlowEdge() : rev(NULL) {}
	FlowEdge(size_t _from, size_t _to, size_t _id, int _flow, int _capacity, FlowEdge* _rev) : Edge(_from, _to, _id) {
		flow = _flow;
	 	capacity = _capacity;
		rev = _rev;
	}
	void setRev(FlowEdge* _rev);
	void push(int x);
	int getG();
	int getFlow();
	int getCapacity();
	FlowEdge* getRev();
	void setFlow(int f);
};

class MKMVertex : public Vertex {
private:
    long long inMaxFlow, outMaxFlow;
    bool deleted;
public:
    MKMVertex(){reset();}
    explicit MKMVertex(Vertex* v) {
        id = v->getId();
        reset();
    }
    void addInFlow(int df);
    void addOutFlow(int df);
    long long getPotential();
    void reset();
    void setDeleted(bool b);
    bool isDeleted();
};

class MKMEdge : public FlowEdge {
private:
    FlowEdge* parent;
    list<Edge*>::iterator inIter, outIter;
public:
    explicit MKMEdge(FlowEdge* edge) {
        parent = edge;
        from = edge->getFrom();
        to = edge->getTo();
        id = edge->getId();
        flow = edge->getFlow();
        capacity = edge->getCapacity();
    }
    long long pushFlow(MKMVertex* from, MKMVertex* to);
    FlowEdge* getParent();
    void setInIter(list<Edge*>::iterator i);
    void setOutIter(list<Edge*>::iterator i);
    list<Edge*>::iterator getInIter();
    list<Edge*>::iterator getOutIter();
};

class Graph {
friend class MKMFinder;
protected:
	size_t EdgeCount;
	vector<list<Edge*> > E;
	vector<Vertex*> V;

public:
	Graph() {}
	Graph(vector<list<Edge*> >& FlowE, const vector<Vertex*>& FlowV, int minG);
	size_t getSize();
};

class FlowGraph : public Graph {
private:
    void push(FlowEdge* e);
	void relabel(size_t a);
	void initPreflow(size_t S);

	void discharge(int u, size_t S, size_t T, queue<size_t>& q, vector<size_t>& used);
public:
	FlowGraph() {}

	int getMaxCap();
	void readGraph(istream& in);
	long long pushFlow(size_t S, long long flow, int minG, size_t T);
	long long Dinic(size_t S, size_t T);
	void printFlow();

	long long RelabelToFront(size_t S, size_t T);
};

class MKMFinder : public Graph {
private:
    vector<list<Edge*> > inE;
    vector<long long> potential;
    FlowGraph* parent;
    size_t S, T;

    void removeEdge(MKMEdge* edge, size_t from, size_t to);
	void pushFlow(vector<list<Edge*> >& edges, size_t S);
	void bfs(FlowGraph* graph, size_t S);
	void buildNetwork(FlowGraph* graph);

	void removeEdges(stack<size_t>& vertices, list<Edge*>& edges);
	void removeVertices();
	long long findMinPotential(MKMVertex*& minVertex);
	void pushFlowToNetwork();
public:
	MKMFinder() {}
	MKMFinder(FlowGraph* graph, size_t S, size_t T) : S(S), T(T) {
        for (size_t i = 0; i < graph->getSize(); ++i) {
            V.push_back(new MKMVertex(graph->V[i]));
        }
        E.assign(graph->getSize(), list<Edge*>());
        inE.assign(graph->getSize(), list<Edge*>());
        bfs(graph, S);
        buildNetwork(graph);
	}
	~MKMFinder() {
        for (size_t i = 0; i < getSize(); ++i) {
            delete V[i];
        }
        for (size_t i = 0; i < getSize(); ++i) {
            for(list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
                delete (*j);
            }
        }
	}

    long long solve();
};

size_t Vertex::getDist() {
	return dist;
}

list<Edge*>::iterator Vertex::getHead() {
	return edgesHead;
}

void Vertex::setDist(size_t x) {
	dist = x;
}

void Vertex::increaseHead() {
	edgesHead++;
}

void Vertex::initHead(list<Edge*>::iterator i) {
	edgesHead = i;
}

long long Vertex::getExcess() {
	return excess;
}

int Vertex::getLabel() {
	return label;
}

void Vertex::setExcess(long long x) {
	excess = x;
}

void Vertex::setLabel(int x) {
	label = x;
}

void Vertex::increaseExcess(long long x) {
	excess += x;
}

void Vertex::setId(size_t x) {
    id = x;
}

size_t Vertex::getId() {
    return id;
}

void Edge::setUsed(bool u) {
    used = u;
}
bool Edge::getUsed() {
    return used;
}

size_t Edge::getTo() {
	return to;
}

size_t Edge::getId() {
	return id;
}

size_t Edge::getFrom() {
	return from;
}

int FlowEdge::getFlow() {
	return flow;
}

void FlowEdge::setRev(FlowEdge* _rev) {
	rev = _rev;
}

void FlowEdge::push(int x) {
	this->flow += x;
	this->rev->flow -= x;
}

int FlowEdge::getG() {
	return capacity - flow;
}

int FlowEdge::getCapacity() {
    return capacity;
}
FlowEdge* FlowEdge::getRev() {
    return rev;
}
void FlowEdge::setFlow(int f) {
    flow = f;
}

void MKMVertex::addInFlow(int df) {
    inMaxFlow += df;
}
void MKMVertex::addOutFlow(int df) {
    outMaxFlow += df;
}
long long MKMVertex::getPotential() {
    return min(inMaxFlow, outMaxFlow);
}
void MKMVertex::reset() {
    inMaxFlow = outMaxFlow = excess = 0;
    deleted = false;
}
void MKMVertex::setDeleted(bool b) {
    deleted = b;
}
bool MKMVertex::isDeleted() {
    return deleted;
}

long long MKMEdge::pushFlow(MKMVertex* from, MKMVertex* to) {
    long long f = min(from->getExcess(), min(to->getPotential(), 0ll + capacity - flow));
    if (from->getId() == getFrom()) {
        from->addOutFlow(-f);
        to->addInFlow(-f);
    } else {
        to->addOutFlow(-f);
        from->addInFlow(-f);
    }
    flow += f;
    to->increaseExcess(f);
    from->increaseExcess(-f);
    return f;
}
FlowEdge* MKMEdge::getParent() {
    return parent;
}
void MKMEdge::setInIter(list<Edge*>::iterator i) {
    inIter = i;
}
void MKMEdge::setOutIter(list<Edge*>::iterator i) {
    outIter = i;
}
list<Edge*>::iterator MKMEdge::getInIter() {
    return inIter;
}
list<Edge*>::iterator MKMEdge::getOutIter() {
    return outIter;
}

size_t Graph::getSize() {
	return E.size();
}

void MKMFinder::bfs(FlowGraph* graph, size_t S) {
	size_t n = graph->E.size();
    for (size_t i = 0; i < n; ++i) {
        V[i]->setDist(INF);
    }
    V[S]->setDist(0);
    queue<size_t> q;
    q.push(S);

    while (!q.empty()) {
        size_t cur = q.front();
        q.pop();

        for (list<Edge*>::iterator i = graph->E[cur].begin(); i != graph->E[cur].end(); ++i) {
            if ((static_cast<FlowEdge*>(*i))->getCapacity() - (static_cast<FlowEdge*>(*i))->getFlow() <= 0)
                continue;
            size_t y = (*i)->getTo();
            if (V[y]->getDist() == INF) {
                V[y]->setDist(V[cur]->getDist() + 1);
                (*i)->setUsed(true);
                q.push(y);
            }
        }
    }
}

void MKMFinder::buildNetwork(FlowGraph* graph) {
    for (size_t i = 0; i < graph->getSize(); ++i) {
        for(list<Edge*>::iterator j = graph->E[i].begin(); j != graph->E[i].end(); ++j) {
            if (V[(*j)->getFrom()]->getDist() != V[(*j)->getTo()]->getDist() - 1)
                continue;
            E[i].push_back(static_cast<Edge*>(new MKMEdge(static_cast<FlowEdge*>(*j))));
            inE[(*j)->getTo()].push_back(E[i].back());
            (static_cast<MKMEdge*>(E[i].back()))->setInIter(--(inE[(*j)->getTo()].end()));
            (static_cast<MKMEdge*>(E[i].back()))->setOutIter((--E[i].end()));
        }
    }
}

Graph::Graph(vector<list<Edge*> >& FlowE, const vector<Vertex*>& FlowV, int minG) {
	E.resize(FlowE.size());
	EdgeCount = 0;
	for (size_t i = 0; i < FlowE.size(); ++i) {
		E[i].resize(0);
		for (list<Edge*>::iterator j = FlowE[i].begin(); j != FlowE[i].end(); ++j) {
			if ((static_cast<FlowEdge*>(*j))->getG() >= minG) {
				size_t y = (*j)->getTo();
				E[i].push_back(new Edge(i, y, (*j)->getId()));
				EdgeCount++;
			}
		}
	}

	V = FlowV;
}

int FlowGraph::getMaxCap() {
	int ans = 0;
	for (size_t i = 0; i < E.size(); ++i) {
		for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
			ans = max(ans, (static_cast<FlowEdge*>(*j))->getFlow() + (static_cast<FlowEdge*>(*j))->getG());
		}
	}
	return ans;
}

void FlowGraph::readGraph(istream& in) {
	int n, m;
	in >> n >> m;
	EdgeCount = m;
	E.resize(n);
	V.resize(n);
	for (size_t i = 0; i < n; ++i) {
		V[i] = new Vertex();
		V[i]->setId(i);
	}

	for (size_t i = 0; i < m; ++i) {
		size_t a, b, c;
		in >> a >> b >> c;
		a--; b--;


		E[a].push_back(new FlowEdge(a, b, i + 1, 0, c, NULL));
		E[b].push_back(new FlowEdge(b, a, 0, 0, 0, static_cast<FlowEdge*>(E[a].back())));
		(static_cast<FlowEdge*>(E[a].back()))->setRev(static_cast<FlowEdge*>(E[b].back()));
	}
}

long long FlowGraph::pushFlow(size_t x, long long flow, int minG, size_t T) {
	if (x == T) {
		return flow;
	}
	long long oldflow = flow;

	for (list<Edge*>::iterator i = V[x]->getHead(); i != E[x].end(); ++i) {
		size_t y = (*i)->getTo();
		if (V[y]->getDist() == V[x]->getDist() + 1 && (static_cast<FlowEdge*>(*i))->getG() >= minG && flow > 0) {
			long long cur = pushFlow(y, min(flow, static_cast<long long>((static_cast<FlowEdge*>(*i))->getG())), minG, T);
			flow -= cur;
			if (flow != 0) {
				V[x]->increaseHead();
			}
			(static_cast<FlowEdge*>(*i))->push(cur);
		}
	}

	return (oldflow - flow);
}

long long FlowGraph::Dinic(size_t S, size_t T) {
	long long flow = 0;
    while (true) {
        MKMFinder net(this, S, T);
        long long cur = net.solve();
        if (cur == 0) {
            break;
        }
        flow += cur;
    }
	return flow;
}


void FlowGraph::push(FlowEdge* e) {
 	size_t a = e->getFrom();
 	size_t b = e->getTo();
 	long long d = min(V[a]->getExcess(), static_cast<long long>(e->getG()));
 	e->push(d);
 	V[a]->increaseExcess(-d);
 	V[b]->increaseExcess(d);
}

void FlowGraph::relabel(size_t a) {
	int ans = INF;
	for (list<Edge*>::iterator i = E[a].begin(); i != E[a].end(); ++i) {
		if ((static_cast<FlowEdge*>(*i))->getG() > 0) {
			ans = min(ans, V[ (*i)->getTo() ]->getLabel() + 1);
		}
	}
	V[a]->setLabel(ans);
}

void FlowGraph::initPreflow(size_t S) {
	for (size_t i = 0; i < V.size(); ++i) {
		V[i]->setExcess(0);
		V[i]->setLabel(0);
	}
	for (list<Edge*>::iterator i = E[S].begin(); i != E[S].end(); ++i) {
		if ((*i)->getTo() == S) {
			continue;
		}
		V[ (*i)->getTo() ]->increaseExcess( (static_cast<FlowEdge*>(*i))->getG() );
		V[S]->increaseExcess(-(static_cast<FlowEdge*>(*i))->getG());
		(static_cast<FlowEdge*>(*i))->push((static_cast<FlowEdge*>(*i))->getG());
	}
	V[S]->setLabel(V.size());
}

void FlowGraph::discharge(int u, size_t S, size_t T, queue<size_t>& q, vector<size_t>& used) {
	while (V[u]->getExcess() > 0) {
		if (V[u]->getHead() == E[u].end()) {
			relabel(u);
			V[u]->initHead(E[u].begin());
		}
		else {
			FlowEdge* e = static_cast<FlowEdge*>(*V[u]->getHead());
			if (e->getG() > 0 && V[u]->getLabel() == V[e->getTo()]->getLabel() + 1) {
				push(e);
				if (used[e->getTo()] == 0 && V[e->getTo()]->getExcess() > 0 && e->getTo() != S && e->getTo() != T) {
					q.push(e->getTo());
					used[e->getTo()] = 1;
				}
			}
			else
				V[u]->increaseHead();
		}
	}
	used[u] = 0;
}

long long FlowGraph::RelabelToFront(size_t S, size_t T) {
	initPreflow(S);
	vector<size_t> used(V.size());
	queue<size_t> q;
	for (size_t i = 0; i < V.size(); ++i) {
		V[i]->initHead(E[i].begin());
		if (i != S && i != T && V[i]->getExcess() > 0) {
			q.push(i);
			used[i] = 1;
		}
	}


	while (q.size() > 0) {
		size_t cur = q.front();
		q.pop();
		discharge(cur, S, T, q, used);
	}

	return V[T]->getExcess();
}

void FlowGraph::printFlow() {
	vector<int> ans;
	ans.resize(EdgeCount);
	for (size_t i = 0; i < E.size(); ++i) {
		for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
			if ((*j)->getId() != 0) {
				size_t id = (*j)->getId();
				ans[id - 1] = (static_cast<FlowEdge*>(*j))->getFlow();
			}
		}
	}

	for (size_t i = 0; i < ans.size(); ++i) {
		cout << ans[i] << endl;
	}
}

void MKMFinder::removeEdge(MKMEdge* edge, size_t from, size_t to) {
    long long df = edge->getCapacity() - edge->getFlow();
    inE[to].erase(edge->getInIter());
    E[from].erase(edge->getOutIter());
    (static_cast<MKMVertex*>(V[from]))->addOutFlow(-df);
    (static_cast<MKMVertex*>(V[to]))->addInFlow(-df);
    edge->getParent()->setFlow(edge->getFlow());
    edge->getParent()->getRev()->setFlow(-edge->getFlow());
    delete edge;
}

void MKMFinder::pushFlow(vector<list<Edge*> >& edges, size_t S) {
    queue<size_t> q;
    q.push(S);
    while(!q.empty()) {
        size_t currNum = q.front();
        q.pop();
        MKMVertex* current = (static_cast<MKMVertex*>(V[currNum]));
        for (list<Edge*>::iterator i = edges[currNum].begin(); i != edges[currNum].end(); ++i) {
            if (current->getExcess() == 0)
                break;
            size_t to = (*i)->getTo();
            if (to == currNum)
                to = (*i)->getFrom();
            (static_cast<MKMEdge*>(*i))->pushFlow(static_cast<MKMVertex*>(V[currNum]), static_cast<MKMVertex*>(V[to]));
            q.push(to);
        }
    }
}

void MKMFinder::removeEdges(stack<size_t>& vertices, list<Edge*>& edges) {
    while (!edges.empty()) {
        MKMEdge* e = static_cast<MKMEdge*>(*edges.begin());
        size_t from = e->getFrom();
        size_t to = e->getTo();
        removeEdge(e, from, to);
        if ((static_cast<MKMVertex*>(V[from]))->getPotential() == 0 && !(static_cast<MKMVertex*>(V[from]))->isDeleted())
            vertices.push(from);
        if ((static_cast<MKMVertex*>(V[to]))->getPotential() == 0 && !(static_cast<MKMVertex*>(V[to]))->isDeleted())
            vertices.push(to);
    }
}

void MKMFinder::removeVertices() {
    stack<size_t> vertices;
    for (size_t i = 0; i < getSize(); ++i) {
        if ((static_cast<MKMVertex*>(V[i]))->getPotential() == 0 && !(static_cast<MKMVertex*>(V[i]))->isDeleted()) {
            vertices.push(i);
        }
    }
    while (!vertices.empty()) {
        int id = vertices.top();
        vertices.pop();
        (static_cast<MKMVertex*>(V[id]))->setDeleted(true);
        removeEdges(vertices, inE[id]);
        removeEdges(vertices, E[id]);

    }
}

long long MKMFinder::findMinPotential(MKMVertex*& minVertex) {
    long long minPotential = INFLL;
    for (size_t i = 0; i < getSize(); ++i) {
        if (!(static_cast<MKMVertex*>(V[i]))->isDeleted() && minPotential > (static_cast<MKMVertex*>(V[i]))->getPotential()) {
            minVertex = (static_cast<MKMVertex*>(V[i]));
            minPotential = minVertex->getPotential();
        }
    }
    return minPotential;
}

void MKMFinder::pushFlowToNetwork() {
    for (size_t i = 0; i < getSize(); ++i) {
        for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
            (static_cast<MKMEdge*>(*j))->getParent()->setFlow((static_cast<MKMEdge*>(*j))->getFlow());
            (static_cast<MKMEdge*>(*j))->getParent()->getRev()->setFlow(-(static_cast<MKMEdge*>(*j))->getFlow());
        }
    }
}

long long MKMFinder::solve() {
    for (size_t i = 0; i < getSize(); ++i) {
        (static_cast<MKMVertex*>(V[i]))->reset();
    }
    for (size_t i = 0; i < getSize(); ++i) {
        for (list<Edge*>::iterator j = E[i].begin(); j != E[i].end(); ++j) {
            (static_cast<MKMVertex*>(V[(*j)->getFrom()]))->addOutFlow((static_cast<MKMEdge*>(*j))->getCapacity() - (static_cast<MKMEdge*>(*j))->getFlow());
            (static_cast<MKMVertex*>(V[(*j)->getTo()]))->addInFlow((static_cast<MKMEdge*>(*j))->getCapacity() - (static_cast<MKMEdge*>(*j))->getFlow());
        }
    }
    (static_cast<MKMVertex*>(V[S]))->addInFlow(INFLL);
    (static_cast<MKMVertex*>(V[T]))->addOutFlow(INFLL);

    long long totalFlow = 0;
    while (true) {
        removeVertices();
        MKMVertex* minVertex;
        long long minPotential = findMinPotential(minVertex);
        if (minPotential == INFLL)
            break;

        minVertex->increaseExcess(minPotential);
        pushFlow(E, minVertex->getId());
        minVertex->increaseExcess(minPotential);
        pushFlow(inE, minVertex->getId());
        totalFlow += minPotential;
    }
    pushFlowToNetwork();
    return totalFlow;
}


int main() {
    //freopen("input.txt", "rt", stdin);
	FlowGraph G;
	G.readGraph(cin);
	cout << G.Dinic(0, G.getSize() - 1) << endl;
	G.printFlow();
	return 0;
}
