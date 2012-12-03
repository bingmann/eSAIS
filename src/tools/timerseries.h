/******************************************************************************
 * src/tools/timerseries.h
 *
 * Simple class to take multiple timer measurements and print out a summary at
 * the end. Also replaceable with a fake class that does nothing.
 *
 ******************************************************************************
 * Copyright (C) 2012 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

class TimerSeries
{
private:
    std::vector<double>		m_times;
    std::vector<std::string>	m_names;

public:

    // take new time record at this position.
    TimerSeries& record(const std::string& name)
    {
	m_times.push_back( omp_get_wtime() );
	m_names.push_back( name );
	return *this;
    }
    
    // take new time record at this position.
    TimerSeries& record() 
    {
	std::ostringstream os;
	os << "ts" << m_times.size() + 1;
	return record(os.str());
    }

    // print out currently last timer mark
    TimerSeries& printlast()
    {
	assert( m_times.size() >= 2 );
	std::cout << m_names[ m_times.size() - 2 ] << "->" << m_names[ m_times.size() - 1 ] << " = "
		  << (m_times[ m_times.size() - 1 ] - m_times[ m_times.size() - 2 ]) << "\n";
	return *this;
    }

    TimerSeries& printall()
    {
	for (unsigned int i = 1; i < m_times.size(); ++i)
	{
	    std::cout << m_names[i-1] << "->" << m_names[i] << " = "
		      << (m_times[i] - m_times[i-1]) << " - ";
	}
	std::cout << "total = " << m_times.back() - m_times.front() << "\n";
	return *this;
    }
};

// Fake class to replace TimerSeries when doing performance
// measurements
class TimerFake
{
public:
    TimerFake& record()
    {
 	return *this;
    }
    TimerFake& record(const std::string& /* name */)
    {
 	return *this;
    }
    TimerFake& printlast()
    {
 	return *this;
    }
    TimerFake& printall()
    {
 	return *this;
    }
};
