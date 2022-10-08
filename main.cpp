#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <sstream>

using namespace std;

#include "show.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//#include <CGAL/Triangulation_3.h>
//#include <CGAL/draw_triangulation_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

typedef Kernel::Point_2 Point2;

typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_set_3<Point> Point_set;

typedef CGAL::Triangulation_2<Kernel>                            Triangulation;
typedef Triangulation::Point                                T_Point;

typedef Triangulation::Vertex_handle 	Vertex_handle;
typedef Triangulation::Face_handle 	Face_handle;

typedef Triangulation::Finite_faces_iterator 	Face_iterator;
typedef Triangulation::Finite_vertices_iterator 	Vertex_iterator;


long double global_z_min, global_z_max;


int filter(string s, vector<string> *elements)
{
    int indicator = 0;
    string element;
    bool in_flag = false, added = false;

    for (int i = 0; i < s.size(); i++)
    {
        if (s[i] != ' ' && s[i] != '\t')
        {
            added = false;
            in_flag =true;
            element += s[i];
            
        }
        else
        {
            in_flag = false;
        }
        
        if (!in_flag && !added)
        {
            added = true;
            elements->push_back(element);
            indicator ++;
            element.clear();
        }
    }
    elements->push_back(element);
    indicator ++;
    return indicator;
    
}

std::string extractLastNChars(std::string const &str, int n)
{
    if (str.size() < n) {
        return str;
    }
 
    return str.substr(str.size() - n);
}

bool read_file(string file_name, Point_set *point_set, vector<Point> *points)
{
    vector<string> elements;
    cout << "Reading file " << file_name << endl;
    if (extractLastNChars(file_name, 4) == ".obj")
    {
        //system("f3d box.obj --up +Z");
        fstream file;
        file.open(file_name, ios::in);

        string inp;
        char del = '\r';

        //First line (info)
        getline(file, inp, del);
        cout << inp;
        //Second line (Date and Time of file creation)
        getline(file, inp, del);
        cout << inp;

        cout << "\n=======================================" << endl;

        stringstream stream;
        long double ac_x, ac_y, ac_z;
        //vector<Point> points;

        int number_of_points=0;
        
        bool first_vertex = true;
        while (getline(file, inp, del))
        {
            elements.clear();
            filter(inp,&elements);      //filter function is not working properly
            //check if readed line is vertex (Starts with v)
            if (elements[0] == "\nv")
            {
                number_of_points ++;
                stream << elements[1];
                stream >> ac_x;
                stream.clear();
                stream << elements[2];
                stream >> ac_y;
                stream.clear();
                stream << elements[3];
                stream >> ac_z;
                stream.clear();
                points->push_back(Point(ac_x, ac_y, ac_z));
                point_set->insert(Point(ac_x, ac_y, ac_z));
                if(first_vertex)
                {
                    global_z_min = global_z_max = ac_z;
                    first_vertex = false;
                }
                else
                {
                    if(ac_z < global_z_min)
                        global_z_min = ac_z;
                    if(ac_z > global_z_max)
                        global_z_max = ac_z;
                }
                
                //cout << "\np" << number_of_points << " = " << ac_x << "\t" << ac_y << "\t" << ac_z << endl;
            }
            
        }
        return true;
    }
    else
    {
        cout << "\nUnsupported file format" << endl;
        return false;
    }
}

void print_point_set (const Point_set& point_set)
{
  cout << "Content of point set:" << endl;
  for (Point_set::const_iterator it = point_set.begin(); it != point_set.end(); ++it)
    cout << "* Point " << *it << " : " << point_set.point(*it) << endl;
}


vector<Point2> xy_projection(const vector<Point> point_vector)
{
    vector<Point2> projection;
    for (int i=0; i<point_vector.size(); i++)
    {
        projection.push_back(Point2(point_vector[i].x() ,point_vector[i].y()));
    }
    return projection;

}


Triangulation triangulation(const vector<Point2> ps)
{
    Triangulation t;
    for (int i=0; i<ps.size(); i++)
    {
      t.insert(T_Point(ps[i].x(),ps[i].y()));
    }

    return t;
}

double distance2(const Point2 p,const Point2 q) {
    return ((p.x() - q.x()) * (p.x() - q.x())) + ((p.y() - q.y()) * (p.y() - q.y()));
}

vector<long double> min_max_z(const Point_set point_set)
{
    bool first_vertex = true;
    long double min,max;
    vector<long double> res;
    for (Point_set::const_iterator it = point_set.begin(); it != point_set.end(); ++it)
    {
        long double ac_z = point_set.point(*it).z();
        if(first_vertex)
        {
            min = max = ac_z;
            first_vertex = false;
        }
        else
        {
            if(ac_z < min)
                min = ac_z;
            if(ac_z > max)
                max = ac_z;
        }
    }
    res.push_back(min);
    res.push_back(max);
    return res;
}

vector<Point_set> slice(const Point_set point_set, int epochs)
{
    vector<Point_set> slices;
    Point_set ac;
    long double k = (global_z_max - global_z_min)/epochs;

    for (int i=0; i < epochs; ++i)
        slices.push_back(ac);


    for (Point_set::const_iterator it = point_set.begin(); it != point_set.end(); ++it)
    {
        for (int i=1; i <= epochs; ++i)
        {
            if(point_set.point(*it).z() > (i - 1)*k && point_set.point(*it).z() <= i*k)
                slices[i-1].insert(point_set.point(*it));
        }            
    }
    return slices;

}

void draw_slices(vector<Point_set> slices, float gap)
{
    Point_set res;
    for (int i = 0; i < slices.size(); ++i)
        for (Point_set::const_iterator it = slices[i].begin(); it != slices[i].end(); ++it)
            res.insert(Point(slices[i].point(*it). x(),slices[i].point(*it).y(), slices[i].point(*it).z() + (i+1)*gap));
    
    CGAL::draw(res);
    
}

vector<Point> point_set_to_vector(Point_set point_set)
{
    vector<Point> points;
    for (Point_set::const_iterator it = point_set.begin(); it != point_set.end(); ++it)
        points.push_back(point_set.point(*it));
    return points;

}

int main(int argc, char** argv)
{
    int number_of_slices = 8;


    if (argc == 1)
    {
        cout << "Enter file path as an argument!" << endl;
        return 1;
    }

    Point_set points;
    vector<Point> point_vector;
    read_file(argv[1], &points, &point_vector);

    //cout << "\nVERTEX : " << point_vector[1].x() << " " << point_vector[1].y() << " " << point_vector[1].z() << endl ; 

    //print_point_set(points);


    CGAL::draw(points);

    
    vector<Point_set> sls = slice(points, number_of_slices);

    draw_slices(sls, 1);


    // calculate area

    Face_handle f = Face_handle();
    //Vertex_handle v = Vertex_handle();


    Face_iterator fit;
    Vertex_iterator vit;
    Triangulation t;
    long double area;
    long double slice_height = (global_z_max - global_z_min) / number_of_slices;
    long double volume = 0;
    

    cout << "Calculating Volume of " << argv[1] << " ..." << endl;
    for (int i=0; i<sls.size(); i++)
    {
        t = triangulation(xy_projection(point_set_to_vector(sls[i])));
        //CGAL::draw(t);
        for (fit = t.finite_faces_begin(); fit != t.finite_faces_end(); fit++)
        {
            // area of each triangle
            area =  CGAL::area(fit->vertex(0)->point(),fit->vertex(1)->point(),fit->vertex(2)->point());

            //area * slice height
            volume += area * slice_height;
        }
    }

    cout << "Volume = " << volume << endl;




}
