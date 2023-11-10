#include <iostream>
#include <stack>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include <cmath>
#include <map>
#include <variant>

class Array {
private:
    int* current_array_;
    int size_;

    static int generate_random_number() {
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        const int random_number = std::rand() % 6 + 5;
        return random_number;
    }

public:
    explicit Array(const int& size) :size_(size)
    {
        current_array_ = new int[size_];
        for (int i = 0; i < size; i++)
        {
            current_array_[i] = generate_random_number();
        }
    }
    Array() = default;

    Array(const Array& other_array) : size_(other_array.size_) {
        current_array_ = new int[size_];
        for (int i = 0; i < size_; i++)
        {
            current_array_[i] = other_array.current_array_[i];
        }
    }

    void ChangeSize(const int& size) {
        const int* tmp_array = current_array_;
        current_array_ = new int[size];

        for (int i = 0; i < size; i++)
        {
            current_array_[i] = tmp_array[i];
        }

        delete[] tmp_array;
        size_ = size;
    }

    void print_array() const
    {
        if (size_ <= 0) return;
        for (size_t i = 0; i < size_; ++i) {
            std::cout << current_array_[i] << " ";
        }
    }
    void sort() const
    {
        if (size_ <= 0) return;
        for (int i = 1; i < size_; i++)
        {
            int j = i - 1;
            while (j >= 0 && current_array_[j] > current_array_[j + 1])
            {
                std::swap(current_array_[j], current_array_[j + 1]);
                j--;
            }
        }
    }

    int find_max() const
    {
        if (size_ <= 0) return -1;
        int max = current_array_[0];
        for (size_t i = 1; i < size_; ++i) {
            if (current_array_[i] > max) {
                max = current_array_[i];
            }
        }
        return max;
    }

    int find_min() const
    {
        if (size_ <= 0) return -1;
        int min = current_array_[0];
        for (size_t i = 1; i < size_; ++i) {
            if (current_array_[i] < min) {
                min = current_array_[i];
            }
        }
        return min;
    }

    int& operator[](const unsigned int index) const
    {
        if (index > size_)
            throw std::exception("Index of Range");
        return current_array_[index];
    }
    int operator()(int value) {
        for (size_t i = 0; i < size_; ++i)
        {
            current_array_[i] + value;
        }
    }

    explicit operator int() const
    {
        int sum = 0;
        for (size_t i = 0; i < size_; ++i)
        {
            sum += current_array_[i];
        }
        return sum;
    }
    operator char* () const
    {
        char* str = new char[size_ + 1];

        for (size_t i = 0; i < size_; ++i)
        {
            str[i] = static_cast<char>(current_array_[i]);
        }
        str[size_] = '\0';

        return str;
    }

    ~Array() {
        delete[] current_array_;
    }
};


template<typename T>
class Queue_
{
public:
    Queue_() :size_(0)
    {
        current_array_ = new T[size_];
        for (int i = 0; i <= size_; i++)
        {
            current_array_[i] = 0;
        }
    }
    void push(T value)
    {
        size_++;
        T* tmp = new T[size_];

        for (int i = 0; i < size_ - 1; i++)
            tmp[i] = current_array_[i];

        tmp[size_ - 1] = value;

        for(size_t i = 0; i< size_; ++i)
        {
            current_array_[i] = tmp[i];
        }

        delete[] tmp;
    }
    T pop()
    {
        T res = current_array_[0];

        size_;
        T* tmp = new T[size_];

        for (int i = 0; i < size_ - 1; i++)
            tmp[i] = current_array_[i];

        for (size_t i = 0; i < size_ + 1; ++i)
        {
            current_array_[i] = tmp[i];
        }

        size_--;

        return res;
    }

    void print() const
    {
	    for(size_t i = 0; i< size_; ++i)
	    {
            std::cout << current_array_[i] << " ";
	    }
    }
private:
    T* current_array_;
    unsigned int size_;
};



int main()
{
    Queue_<int> que;
    que.push(5);
    que.push(2);
    que.print();
    std::cout << "\n" <<  que.pop();
    que.print();
}

template<typename T>
class Var
{
public:
    Var(const T variable) : variable_(variable) {}

    Var operator +(const Var& other)
    {
        if(typeid(other.variable_) != typeid(std::string) && typeid(other.variable_) != typeid(char))
        {
            Var result;
            result.variable_ = variable_ + other.variable_;
            return result;
        }

        if(typeid(other.variable_) == typeid(std::string) || typeid(other.variable_) == typeid(char))
        {
            Var result;
            result.variable_ = variable_ + static_cast<int>(other.variable_);
            return result;
        }
    }
    Var operator +=(const Var& other)
    {
        if (typeid(other.variable_) != typeid(std::string) && typeid(other.variable_) != typeid(char))
        {
            variable_ += other.variable_;
            return *this;
        }
        if(typeid(other.variable_) == typeid(std::string) || typeid(other.variable_) == typeid(char))
        {
            variable_ += static_cast<int>(other.variable_);
            return *this;
        }
    }

    Var operator -(const Var& other)
    {
        if (typeid(other.variable_) != typeid(std::string) && typeid(other.variable_) != typeid(char))
        {
            Var result;
            result.variable_ = variable_ - other.variable_;
            return result;
        }

        if (typeid(other.variable_) == typeid(std::string) || typeid(other.variable_) == typeid(char))
        {
            Var result;
            result.variable_ = variable_ - static_cast<int>(other.variable_);
            return result;
        }
    }
    Var operator -=(const Var& other)
    {
        if (typeid(other.variable_) != typeid(std::string) && typeid(other.variable_) != typeid(char))
        {
            variable_ -= other.variable_;
            return *this;
        }
        if (typeid(other.variable_) == typeid(std::string) || typeid(other.variable_) == typeid(char))
        {
            variable_ -= static_cast<int>(other.variable_);
            return *this;
        }
    }

    Var operator *(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            Var result;
            const std::string& str1 = variable_;
            const std::string& str2 = other.variable_;

            for (char c : str1) {
                if (str2.find(c) != std::string::npos) {
                    result += c;
                }
            }

            return result;
        }
        Var result;
        result.variable_ = variable_ * other.variable_;

        return result;
    }
    Var operator *=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            Var result;
            const std::string& str1 = variable_;
            const std::string& str2 = other.variable_;

            for (char c : str1) {
                if (str2.find(c) != std::string::npos) {
                    result += c;
                }
            }

            return result;
        }
        variable_ *= other.variable_;
        return *this;
    }

    Var operator /(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            Var result;
            const std::string& str1 = variable_;
            const std::string& str2 = other.variable_;

            for (char c : str1) {
                if (str2.find(c) == std::string::npos) {
                    result += c;
                }
            }

            return result;
        }

        if(other.variable_ != 0)
        {
            Var result;
            result.variable_ = variable_ / other.variable_;
            return result;
        }else
        {
            throw std::exception("Divide by zero");
        }
    }
    Var operator /=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            Var result;
            const std::string& str1 = variable_;
            const std::string& str2 = other.variable_;

            for (char c : str1) {
                if (str2.find(c) == std::string::npos) {
                    result += c;
                }
            }

            return result;
        }

        if (other.variable_ != 0)
        {
            variable_ /= other.variable_;
            return *this;
        }
        else throw std::exception("Divide by zero");
    }

    bool operator <(const Var& other)
    {
	    if(typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
	    {
            return variable_.length() < other.variable_.length();
	    }else
	    {
            return variable_ < other.variable_;
	    }
    }
    bool operator <=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() <= other.variable_.length();
        }
        else
        {
            return variable_ <= other.variable_;
        }
    }

    bool operator >(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() > other.variable_.length();
        }
        else
        {
            return variable_ > other.variable_;
        }
    }
    bool operator >=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() >= other.variable_.length();
        }
        else
        {
            return variable_ >= other.variable_;
        }
    }

    bool operator ==(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() == other.variable_.length();
        }
        else
        {
            return variable_ == other.variable_;
        }
    }

    bool operator !=(const Var& other)
    {
        if (typeid(T) == typeid(std::string) || typeid(T) == typeid(char))
        {
            return variable_.length() != other.variable_.length();
        }
        else
        {
            return variable_ != other.variable_;
        }
    }

private:
    T variable_;
};

template<typename T>
class Arr
{
private:
    T* current_array_;
    unsigned int size_;

    static T generate_random_number() {
        std::srand(static_cast<T>(std::time(nullptr)));
        const T random_number = std::rand() % 100 + 1;
        return random_number;
    }

public:
    explicit Arr(const unsigned int& size) :size_(size)
    {
        current_array_ = new T[size_];
        for (size_t i = 0; i < size; ++i)
        {
            current_array_[i] = generate_random_number();
        }
    }
    Arr() = default;

    T get_value(const unsigned int& index) const
    {
	    if(index > size_)
	    {
            return throw std::out_of_range("Bad index");
	    }
        return *current_array_[index];
    }

    unsigned int get_size() const { return  size_; }

    void push_value(const unsigned int& index,const T& value)
    {
	    if(index > size_)
	    {
            return throw std::out_of_range("Bad index");
	    }

        // Надо организовать какую-то проверку указателя массива T и типа value

        T* tmp_array_ = new T[size_ + 1];

        for(size_t i = 0; i< size_; ++i)
        {
            if(i == index)
            {
                /*
					* class Foo{
					* public:
					*   Foo(const int value) : value(value) {}
					* };
					* ***** Конструктора по умолчанию нет!
					* ***** Если использоловать простое присваиваение *tmp_array_[i] = value;* ,то выпадет исключение т.к будет попытка вызова конструктора по умолчанию, которого нет
					* ***** Код - *new (&tmp_array_[i]) T(value);*  Вызовет конструктор копирования(будь то неявный(дефолтный), || явный(описанный));
				*/
                new (&tmp_array_[i]) T(value);
                //tmp_array_[i] = value;
                continue;
            }
            tmp_array_[i] = current_array_[i];
        }

        delete[] current_array_;
        current_array_ = tmp_array_;
        ++size_;
    }

    void print_array() const
    {
	    for(auto element : current_array_)
	    {
            std::cout << element << " ";
	    }
    }
};

template<typename T>
class Linear
{
private:
    T valueA_;
    T valueX_;
    T valueB_;
    T valueC_;
public:
    Linear(T valueA_, T valueX_, T valueB_)
	:  valueA_(valueA_), valueX_(valueX_), valueB_(valueB_) {}

    Linear(T valueA_, T valueX_, T valueB_, T valueC_)
        : valueA_(valueA_), valueX_(valueX_), valueB_(valueB_), valueC_(valueC_) {}


};

template<typename T>
T get_square(T value_a, T value_b)
{
    T sq;
    sq = -(value_b)/value_a;

    return sq;
}
template<typename T, int>
auto get_square(T value_a, T value_b, T value_c)
{
    T D = std::pow(value_b,2) - 4 * value_a * value_c;

    if(D < 0)
    {
        return nullptr;
    }
    if(D == 0)
    {
        T root = -value_b / (2 * value_a);
        return root;
    }
    if (D > 0)
    {
        T root_one = (-value_b + std::sqrt(D)) / (2 * value_a);
        T root_two = (-value_b - std::sqrt(D)) / (2 * value_a);
        std::pair<T, T> sq = { root_one, root_two };

        return sq;
    }
}

template<typename T>
class Matrix
{
public:
    explicit Matrix(const unsigned int& rows, const unsigned int& columns) : rows_(rows), columns_(columns)
    {
        matrix_ = new T * [rows_];
        for (size_t i = 0; i < rows_; i++)
        {
            matrix_[i] = new T[columns_];
        }
    }
    unsigned int get_row() const
    {
        return rows_;
    }
    unsigned int get_colum() const
    {
        return columns_;
    }

    void fillFromKeyboard()
    {
        for (size_t i = 0; i < rows_; i++)
        {
            for (size_t j = 0; j < columns_; ++j)
            {
                std::cin >> matrix_[i][j];
            }
        }
    }
    void fillRandom()
    {
        for (size_t i = 0; i < rows_; i++)
        {
            for (size_t j = 0; j < columns_; ++j)
            {
                matrix_[i][j] = generate_random_number();
            }
        }
    }

    Matrix operator + (const Matrix& other)
    {
	    if(this.get_row() != other.get_row() || this.get_colum() != other.get_colum())
	    {
            throw std::invalid_argument("The matrix must have the same size for addition.");
	    }
        unsigned int rows = this.get_row();
        unsigned int columns = this.get_colum();
        Matrix<T> result(rows, columns);

        for(size_t i = 0; i < rows; ++i)
        {
	        for(size_t j = 0; j < columns; ++j)
	        {
                result[i][j] = this->matrix_[i][j] + other.matrix_[i][j];
	        }
        }

        return result;
    }
    Matrix operator - (const Matrix& other)
    {
        if (this.get_row() != other.get_row() || this.get_colum() != other.get_colum())
        {
            throw std::invalid_argument("The matrix must have the same size for subtraction.");
        }
        unsigned int rows = this.get_row();
        unsigned int columns = this.get_colum();
        Matrix<T> result(rows, columns);

        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = 0; j < columns; ++j)
            {
                result[i][j] = this.matrix_[i][j] - other.matrix_[i][j];
            }
        }

        return result;
    }
    T* operator[](size_t row_index)
    {
        if (row_index < rows_)
        {
            return matrix_[row_index];
        }
        throw std::out_of_range("Row index out of range");
    }
    Matrix operator * (const Matrix& other)
    {
	    if(this->get_colum() != other.get_colum())
	    {
            throw std::invalid_argument("The number of columns in the matrix must be the same for multiplication");
	    }
	    const unsigned int rows = this->get_row();
	    const unsigned int columns = other.get_colum();
	    const unsigned int common_dimension = this->get_colum();
        unsigned int result_columns;

        if (columns < common_dimension)
        {
            result_columns = common_dimension;
        }
        else
        {
            result_columns = columns;
        }
        Matrix<T> result(rows, result_columns);

        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = 0; j < columns; ++j)
            {
                result[i][j] = T();
                for (size_t k = 0; k < common_dimension; ++k)
                {
                    result[i][j] += this->matrix_[i][k] * other.matrix_[k][j];
                }
            }
        }
        return result;
    }
    
    ~Matrix()
    {
        for (unsigned int i = 0; i < rows_; ++i)
        {
            delete[] matrix_[i];
        }
        delete[] matrix_;
    }
private:
    T** matrix_;
    unsigned int rows_ = 0;
    unsigned int columns_ = 0;

    static T generate_random_number() {
        if(std::is_arithmetic<T>::value)
        {
            std::srand(static_cast<unsigned>(std::time(nullptr)));
            const T random_number = std::rand() % 100 + 1;
            return random_number;
        }else
        {
            std::srand(static_cast<unsigned>(std::time(nullptr)));
            const T random_symol = std::rand() % 205 + 50;
            return random_symol;
        }
    }
};

template<typename T>
class Vector
{
public:
    explicit Vector(const unsigned int& size, const bool fill) : current_array_(nullptr), size_(size)
    {
	    if (fill)
	    {
		    for (size_t i = 0; i < size_; i++)
		    {
			    current_array_[i] = generate_random_value();
		    }
	    }
	    else
	    {
		    for (size_t i = 0; i < size_; i++)
		    {
			    current_array_[i] = 0;
		    }
	    }
    }

    unsigned int get_size() const
    {
        return size_;
    }

    //SetSize(int size, int grow = 1) – установка размера массива
    //(если параметр size больше предыдущего размера массива, то выделяется дополнительный блок памяти, если нет, то «лишние» элементы теряются и память освобождается); параметр
	//grow определяет для какого количества элементов необходимо выделить память, если количество элементов превосходит
	//текущий размер массива.Например, SetSize(5, 5); означает, что при добавлении 6 - го элемента размер массива становится
	//равным 10, при добавлении 11 - го - 15 и т.д.;

	void set_size(const unsigned int& size, const unsigned int& grow = 1)
    {
	    if(size > this->size_)
	    {
            T* tmp_array = new T[size];

            for (size_t i = 0; i < size; ++i)
            {
                if (i > size_) break;
                tmp_array[i] = current_array_[i];
            }

            delete[] current_array_;
            current_array_ = nullptr;

            size_ = size;
            current_array_ = new T[size_];

            for (size_t i = 0; i < size_; ++i)
            {
                current_array_[i] = tmp_array[i];
            }

            delete[] tmp_array;
            tmp_array = nullptr;

	    }else
	    {
            T* tmp_array = new T[size];
            size_ = size;

            for(size_t i = 0; i < size; ++i)
            {
                tmp_array[i] = current_array_[i];
            }

            delete[] current_array_;
            current_array_ = nullptr;

            size_ = size;
            current_array_ = new T[size_];

            for (size_t i = 0; i < size_; ++i)
            {
                current_array_[i] = tmp_array[i];
            }

            delete[] tmp_array;
            tmp_array = nullptr;
	    }
    }

private:
    static T generate_random_value() {
        if (std::is_arithmetic<T>::value)
        {
            std::srand(static_cast<unsigned>(std::time(nullptr)));
            const T random_number = std::rand() % 100 + 1;

            return random_number;
        }

    	std::srand(static_cast<unsigned>(std::time(nullptr)));
    	const T random_symol = std::rand() % 205 + 50;

    	return random_symol;
    }

    T* current_array_;
    unsigned int size_;

};

class Stack_
{
public:
	explicit Stack_(const int& size) : max_size_(size), top_index_(-1) {
        stack_array_ = new char[max_size_];
    }

    bool is_full() const {
        return top_index_ == max_size_ - 1;
    }
    bool is_empty() const {
        return top_index_ == -1;
    }

    int get_size() const
	{
        return top_index_ +1;
	}

    void push(const char c) {
        if (is_full()) {
            throw std::overflow_error("Stack is full");
        }
        stack_array_[++top_index_] = c;
    }
    char pop()
	{
		if(is_empty())
		{
            throw std::overflow_error("Stack is empty");
		}
        return stack_array_[top_index_--];
	}
    void pop_c()
	{
        if (is_empty())
        {
            throw std::overflow_error("Stack is empty");
        }
        stack_array_[top_index_--];
	}

    char pop_non_pop() const
    {
		if(is_empty())
		{
            throw std::overflow_error("Stack is empty");
		}
        return stack_array_[top_index_];
	}

    void clear()
	{
		while (!is_empty())
		{
            pop_c();
		}
	}

    ~Stack_()
    {
        delete[] stack_array_;
    }
private:
    char* stack_array_;
    int max_size_;
    int top_index_;
};

void check_balance(const std::string& str)
{
    unsigned int size = str.size();
    Stack_ stack(size);
    std::stack<char> bracketStack;
    int parenthesis_bracket = 0;
    int square_bracket = 0;
    int curly_brace = 0;

    for (size_t i = 0; i < size; ++i)
    {
        stack.push(str[i]);
    }

    while (!stack.is_empty())
    {
        char symbol = stack.pop();
        if (symbol == '(' || symbol == ')' || symbol == '[' || symbol == ']' || symbol == '{' || symbol == '}')
        {
            bracketStack.push(symbol);
        }
    }

    while (!bracketStack.empty())
    {
        char top = bracketStack.top();
        bracketStack.pop();

        if (top == '(' || top == ')')
        {
            if (top == '(')
                ++parenthesis_bracket;
            else
                --parenthesis_bracket;
        }
        else if (top == '[' || top == ']')
        {
            if (top == '[')
                ++square_bracket;
            else
                --square_bracket;
        }
        else if (top == '{' || top == '}')
        {
            if (top == '{')
                ++curly_brace;
            else
                --curly_brace;
        }

    }

    if (parenthesis_bracket == 0 && square_bracket == 0 && curly_brace == 0)
    {
        std::cout << "Balanced";
    }
    else
    {
        std::cout << "Non balanced";
    }
}

bool isBalanced(std::string input) {
	std:: stack<char> stk;

    for (char c : input) {
        if (c == '(' || c == '{' || c == '[') {
            stk.push(c);
        }
        else if (c == ')' || c == '}' || c == ']') {
            if (stk.empty()) {
                return false;
            }

            char top = stk.top();
            stk.pop();

            if ((c == ')' && top != '(') || (c == '}' && top != '{') || (c == ']' && top != '[')) {
                return false;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PeopleINFO
{
public:
    PeopleINFO(const  std::string& name, const  std::string& lastname, const std::string& patronymic, const uint16_t& age)
        :name_(name), lastname_(lastname), patronymic_(patronymic), age_(age)
	{}

    void set_name(const  std::string& name)
    {
        name_ = name;
    }
    void set_lastname(const  std::string& lastname)
    {
        lastname_ = lastname;
    }
    void set_patronymic(const  std::string& patronymic)
    {
        patronymic_ = patronymic_;
    }

    std::string get_name() const
    {
        return name_;
    }
    std::string get_lastname() const
    {
        return lastname_;
    }
    std::string get_patronymic() const
    {
        return patronymic_;
    }

    uint16_t get_age() const
    {
        return age_;
    }

private:
    std::string name_;
    std::string lastname_;
    std::string patronymic_;
    uint16_t age_;

};
class Passport
{
public:
    Passport(const std::string& name,
        const  std::string& lastname,
        const std::string& patronymic,
        const uint16_t& age,
        const std::string& place_registration = "")
        :people_(name, lastname, patronymic, age), series_pass_(rand() % 9999 + 1000), number_pass_(rand() % 999999 + 100000), place_registration_(place_registration)
	{}

    void change_name(const  std::string& name)
    {
        people_.set_name(name);
    }
    void change_lastname(const  std::string& lastname)
    {
        people_.set_lastname(lastname);
    }
    void change_patronymic(const  std::string& patronymic)
    {
        people_.set_patronymic(patronymic);
    }

    void set_series_pass(const unsigned int series_pass)
    {
        series_pass_ = series_pass;
    }
    void set_number_pass(const unsigned int number_pass)
    {
        number_pass_ = number_pass;
    }
    void change_place_registration(const std::string& place_registration)
    {
        place_registration_ = place_registration_;
    }

    std::string get_name() const
    {
        return people_.get_name();
    }
    std::string get_lastname() const
    {
        return people_.get_lastname();
    }
    std::string get_patronymic() const
    {
        return people_.get_patronymic();
    }
    uint16_t get_age() const
    {
        return people_.get_age();
    }

    unsigned int get_series_pass() const
    {
        return series_pass_;
    }
    unsigned int get_number_pass() const
    {
        return number_pass_;
    }
    std::string get_place_registration() const
    {
        return place_registration_;
    }

    void print_info() const
    {
        std::cout
            << "Name: " << get_name() << std::endl
            << "Lastname: " << get_lastname() << std::endl
            << "Patronymic" << get_patronymic() << std::endl
            << "Age: " << get_age() << std::endl
            << "Series Pass: " << get_series_pass() << std::endl
            << "Number Pass: " << get_number_pass() << std::endl
            << "Place Registration: " << get_place_registration() << std::endl;
    }

protected:
    PeopleINFO people_;
    unsigned int series_pass_;
    unsigned int number_pass_;
    std::string place_registration_;
};
class InternationalPassport : public Passport
{
public:
    InternationalPassport(const std::string& name, const  std::string& lastname, const std::string& patronymic, const uint16_t& age):
	    Passport(name, lastname, patronymic, age)
    {}

    void add_visa(std::string& VISA)
    {
        VISA_.emplace_back(VISA);
    }

    std::vector<std::string> get_visa() const
    {
        return VISA_;
    }
    void print_info() const
    {
        print_info();
        std::cout << "VISA: ";
        for (const auto& i : VISA_)
        {
            std::cout << i << ", ";
        }
    }

private:
    std::vector<std::string> VISA_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Circle
{
public:
    Circle(const unsigned int& diameter) : diameter_(diameter) {}

    unsigned int get_diameter() const { return  diameter_; }
private:
    unsigned int diameter_;
};
class Square
{
public:
	Square(const unsigned int& side) :side_(side) {}

    unsigned int get_side() const { return  side_; }
private:
    unsigned int side_;
};
class CircleInSquare :public Circle, public  Square
{
public:
    CircleInSquare(const unsigned int& diameter, const unsigned int& side) : Circle(diameter), Square(side) {}

    std::string isCorrect()
    {
	    if(get_diameter() > get_side())
	    {
            isCorrect_ = "NO";
            return isCorrect_;
	    }
        isCorrect_ = "YES";
        return isCorrect_;
    }
private:
    std::string isCorrect_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Engine
{
public:
	Engine(const unsigned int& power): power_(power) {}

    unsigned int get_power() const { return  power_; }

private:
    unsigned int power_;
};
class Clutch
{
public:
    Clutch(const unsigned int& ratio) :ratio_(ratio) {}

    unsigned int get_ratio() const { return  ratio_; }
private:
    unsigned int ratio_;
};
class Wheel
{
public:
	Wheel(const unsigned int& pressure) :pressure_(pressure) {}

    unsigned int get_pressure() const { return  pressure_;  }
private:
    unsigned int pressure_;
};
class Car : public Engine, public Clutch, public Wheel
{
public:
	Car(const unsigned int& power, const unsigned int& ratio, const unsigned int& pressure, const std::string& mark) :Engine(power), Clutch(ratio), Wheel(pressure), mark_(mark) {}

    std::string get_mark() const { return mark_; }

    void get_info() const
	{
        std::cout
            << "Car mark: " << get_mark() << std::endl
            << "Engine power: " << get_power() << std::endl
            << "Clutch: " << get_ratio() << std::endl
            << "Wheel: " << get_pressure() << std::endl;
	}
private:
    std::string mark_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AcademicGroup
{
public:
	AcademicGroup(const std::string& name): name_academic_(name) {}

    std::string get_name_academic() const
	{
        return name_academic_;
	}

private:
    std::string name_academic_;
};
class Specialization
{
public:
    Specialization(const std::string& name) : name_specialization_(name) {}

    std::string get_name_secialization() const
    {
        return name_specialization_;
    }
private:
    std::string name_specialization_;
};
class  Discipline
{
public:
    Discipline(const std::string& name) : name_descipline_(name) {}

    std::string get_name_descipline() const
    {
        return name_descipline_;
    }
private:
    std::string name_descipline_;
};
class Auditorium
{
public:
    Auditorium(const unsigned int number) : number_auditorium_(number) {}

    unsigned int get_number_auditorium_() const
    {
        return number_auditorium_;
    }
private:
    unsigned int number_auditorium_;
};
class Lecturer
{
public:
    Lecturer(const std::string& name) : name_lecturer_(name) {}

    std::string get_name_lecturer() const
    {
        return name_lecturer_;
    }
private:
    std::string name_lecturer_;
};
class Organizer : public AcademicGroup, public Specialization, public Discipline, public Auditorium, public  Lecturer
{
public:
    void print_info() const
    {
        std::cout
            << "Academic group: " << get_name_academic() << std::endl
            << "Specialization: " << get_name_secialization() << std::endl
            << "Discipline: " << get_name_descipline() << std::endl
            << "Auditorium: " << get_number_auditorium_() << std::endl
            << "Name Lecturer: " << get_name_lecturer() << std::endl;
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Pet
{
public:
	virtual ~Pet() = default;

	Pet(const std::string pet_name) : pet_name_(pet_name) {}

    virtual void Sound() const = 0;
    virtual void Show() const = 0;
    virtual void Type() const = 0;

protected:
    std::string pet_name_;
};
class Dog : public Pet {
public:
    Dog(const std::string dog_name) : Pet(dog_name), dog_name_(dog_name) {}

    virtual void Sound() const override {
        std::cout << "Gav Gav Gav" << std::endl;
    }
    virtual void Show() const override {
        std::cout << "My name is " << dog_name_ << std::endl;
    }
    virtual void Type() const override {
        std::cout << "I'm a dog" << std::endl;
    }

private:
    std::string dog_name_;
};
class Cat : public Pet {
public:
    Cat(const std::string cat_name) : Pet(cat_name), cat_name_(cat_name) {}

    virtual void Sound() const override {
        std::cout << "Mew mew mew" << std::endl;
    }
    virtual void Show() const override {
        std::cout << "My name is " << cat_name_ << std::endl;
    }
    virtual void Type() const override {
        std::cout << "I'm a cat" << std::endl;
    }

private:
    std::string cat_name_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class LOGText
{
public:
    LOGText(const std::string path) :path_(path) {}

    void setPath(const std::string& path)
	{
        path_ = path;
	}
    std::string getPath() const
    {
        return path_;
    }

    virtual void set_content(std::string content) = 0;

    virtual void display() const = 0;
protected:
    std::string path_;
};
class LOGASKIIConcole : LOGText
{
public:
    LOGASKIIConcole(const std::string path, const std::string& text) : LOGText(path) , text_(text) {}

    virtual  void display() const override
    {
        if (content_.empty()) return;
        for (const char i : content_)
        {
            std::cout << static_cast<int>(i);
        }
    }
    virtual void set_content(std::string content) override
    {
        content_ = content;
    }
private:
    std::string text_;
    std::string content_;

    std::string readFile()
    {
        std::string line;

        std::ifstream file(getPath());
        if (file.is_open())
        {
            while (std::getline(file, line))
            {
                content_ += line;
            }
        }
    }
};
class LOGbinaryConcole : LOGText
{
public:
    LOGbinaryConcole(const std::string path, const std::string& text) : LOGText(path), text_(text) {}

    virtual  void display() const override
    {
        if (content_.empty()) return;
        for (char i : content_)
        {
            std::cout << std::bitset<8>(i) << ' ';
        }
    }
    virtual void set_content(std::string content) override
    {
        content_ = content;
    }

private:
    std::string text_;
    std::string content_;

    std::string readFile()
    {
        std::string line;

        std::ifstream file(getPath());
        if (file.is_open())
        {
            while (std::getline(file, line))
            {
                content_ += line;
            }
        }
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Employer
{
public:
    virtual void print() const = 0;

    virtual ~Employer() = default;
};
class President : public Employer
{
public:
	President(const std::string& rang, const unsigned int& salary) :rang_(rang), salary_(salary) {}
    virtual  void print() const override
    {
        std::cout << "I;m " << rang_ << " My salary: " << salary_ << std::endl;
    }
private:
    unsigned int salary_;
    std::string rang_;
};
class Manager : public Employer
{
public:
	Manager(const std::string& rang, const unsigned int& salary) :rang_(rang), salary_(salary) {}

    virtual  void print() const override
	{
        std::cout << "I;m " << rang_ << " My salary: " << salary_ << std::endl;
	}
private:
    unsigned int salary_;
    std::string rang_;
};
class Worker : public Employer
{
public:
    Worker(const std::string& rang, const unsigned int& salary) :rang_(rang), salary_(salary) {}

    virtual  void print() const override
    {
        std::cout << "I;m " << rang_ << " My salary: " << salary_ << std::endl;
    }
private:
    unsigned int salary_;
    std::string rang_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Shape
{
public:
    virtual double_t getSqare() = 0;
};
class Triangle : public Shape
{
public:
	Triangle(const double_t& height, const double_t& footing) :height_(height), footing_(footing) {}

    virtual double_t getSqare() override
	{
        sqare_ = 0.5 * (footing_ * height_);

        return sqare_;
	}
private:
    double_t height_;
    double_t footing_;
    double_t sqare_;
};
class Circle_ : public  Shape
{
public:
    Circle_(const double_t& radius) :radius_(radius) {}

    virtual double_t getSqare() override
    {
        sqare_ = std::pow(radius_, 2) * 3.14;
        return  sqare_;
    }
private:
    double_t radius_;
    double_t sqare_;
};
class Trapezoid : public Shape
{
public:
    Trapezoid(const double_t& height, const double_t& footing_UP, const double_t& footing_DOWN) :height_(height), footing_UP_(footing_UP), footing_DOWN_(footing_DOWN) {}

    virtual double_t getSqare() override
    {
        sqare_ = ((footing_UP_ + footing_DOWN_) / 2) * height_;

        return  sqare_;
    }
private:
    double_t height_;
    double_t footing_DOWN_;
    double_t footing_UP_;
    double_t sqare_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Equalization
{
public:
    virtual void print_result() = 0;
};
class LinearEqualization : public  Equalization
{
public:
    LinearEqualization(const double_t& value_a, const double_t& value_b) : value_a_(value_a), value_b_(value_b) {}

    virtual void print_result() override
    {
        std::cout << "Result Linear Equalization = " << decide() << std::endl;
    }

private:
    double_t value_a_;
    double_t value_b_;

    double_t decide()
    {
        return -(value_b_) / value_a_;
    }
};
class SquareEqualization : public Equalization
{
public:
    SquareEqualization(const double_t& value_a, const double_t& value_b, const double_t value_c) : value_a_(value_a), value_b_(value_b), value_c_(value_c) {}

    virtual void print_result() override
    {
        std::variant <double_t, std::pair<double_t, double_t>> results = decide();

        if (results.index() == 0)
        {
            double_t root = std::get<double_t>(results);
            if (root > 0)
            {
                std::cout << "Result Square Equalization = " << root;
            }
            else
            {
                std::cout << "Not roots";
            }
        }
        else if (results.index() == 1)
        {
            std::pair<double_t, double_t> roots = std::get<std::pair<double_t, double_t>>(results);
            std::cout << "Root 1 = " << roots.first << " " << "Root 2 = " << roots.second << std::endl;
        }
    }
private:
    double_t value_a_;
    double_t value_b_;
    double_t value_c_;

    std::variant<double_t, std::pair<double_t, double_t>> decide()
    {
        double_t D = std::pow(value_b_, 2) - 4 * value_a_ * value_c_;

        if (D < 0)
        {
            //return nullptr;
        }
        if (D == 0)
        {
            double_t root = -value_b_ / (2 * value_a_);
            return root;
        }
        if (D > 0)
        {
            double_t root_one = (-value_b_ + std::sqrt(D)) / (2 * value_a_);
            double_t root_two = (-value_b_ - std::sqrt(D)) / (2 * value_a_);
            std::pair<double_t, double_t> square = { root_one, root_two };

            return square;
        }
    }
};


class Shape_
{
public:
    virtual void Show() const = 0;
    virtual void Save() const = 0;
    virtual void Load() const = 0;

    virtual ~Shape_() = default;
};
class Square_ : public Shape_
{
public:
	Square_(const double_t& side) :side_(side) {}

    virtual void Show() const override
    {
        std::cout << "sq: " << side_ << ", " << side_ << std::endl;
    }
    virtual void Save() const override
    {
        std::ofstream out;
        out.open("file.txt");
        if(out.is_open())
        {
            out  << "sq: " << side_ << ", " << side_ << std::endl;
        }
        out.flush();
        out.close();
    }
    virtual void Load() const override
    {
        std::string line;
        std::ifstream in;
        in.open("file.txt");
        if(in.is_open())
        {
            while (std::getline(in, line))
            {
                if(line.find("sq: ") != std::string::npos)
					std::cout << line << std::endl;
            }
        }
        in.close();
    }
private:
    double_t side_;
    std::vector<double_t> coordinates_;

    std::vector<double_t> CalclCords()
    {
	    
    }
};
class _Circle_ : public Shape_
{
public:
	_Circle_(const std::double_t radius) : radius_(radius) {}

    virtual void Show() const override
    {
		std::cout << "cir: " << radius_ << ", " << "Center: " << radius_ / 2.0 << std::endl;
    }
    virtual void Save() const override
    {
        std::ofstream out;
        out.open("file.txt");
        if (out.is_open())
        {
            out << "cir: " << radius_ << ", " << "Center: " << radius_ / 2.0 << std::endl;
        }
        out.flush();
        out.close();
    }
    virtual void Load() const override
    {
        std::string line;
        std::ifstream in;
        in.open("file.txt");
        if (in.is_open())
        {
            while (std::getline(in, line))
            {
                if (line.find("cir: ") != std::string::npos)
                    std::cout << line << std::endl;
            }
        }
        in.close();
    }
private:
	double_t radius_;
};
class Ellipse : public Shape_
{
public:
	Ellipse(const double_t radius_a, const double_t radius_b) : radius_a_(radius_a), radius_b_(radius_b) {}

    virtual void Show() const override
    {
        std::cout << "eli: " << "RaA: " << radius_a_ << ", " << "RaB: " << radius_b_ << " Center" << radius_a_/2.0 << std::endl;
    }
    virtual void Save() const override
    {
        std::ofstream out;
        out.open("file.txt");
        if (out.is_open())
        {
            out << "eli: " << "RaA: " << radius_a_ << ", " << "RaB: " << radius_b_ << " Center" << radius_a_ / 2.0 << std::endl;
        }
        out.flush();
        out.close();
    }
    virtual void Load() const override
    {
        std::string line;
        std::ifstream in;
        in.open("file.txt");
        if (in.is_open())
        {
            while (std::getline(in, line))
            {
                if (line.find("eli: ") != std::string::npos)
                    std::cout << line << std::endl;
            }
        }
        in.close();
    }
private:
    double_t radius_a_;
    double_t radius_b_;
};

class Studet
{
public:
    Studet(){}
	Studet(const std::string& name, const std::string& lastname, const std::string& name_group) : name_(name), lastname_(lastname), name_group_(name_group) {}

    std::string get_name() const { return  name_; }
    std::string get_lastname() const { return  lastname_; }
    std::string get_name_group() const { return  name_group_; }

    bool put_estimation(std::string object, uint16_t value)
	{
        if (estimations_.find(object) != estimations_.end())
        {
            std::vector<uint16_t> estimation = estimations_.operator[](object);
            estimation.emplace_back(value);
            estimations_.operator[](object) = estimation;
        }else
        {
            std::map<std::string, std::vector<uint16_t>> my_map;
            std::vector<uint16_t> estimation; estimation.emplace_back(value);

            std::pair<std::string, std::vector<uint16_t>> element_to_insert(object, estimation);
            my_map.insert(element_to_insert);
            estimations_ = my_map;
        }
        return true;
	}

    std::map<std::string, std::vector<uint16_t>> get_estimations() const { return  estimations_; }

protected:
    std::string name_;
    std::string lastname_;
    std::string name_group_;
    std::map<std::string, std::vector<uint16_t>> estimations_;
};
class Group
{
public:
	Group(const std::string& name_group) : name_group_(name_group) {}

    void add_student(Studet studet)
	{
        group_.emplace_back(studet);
	}

    void print() const
	{
        for(size_t i = 0; i < group_.size(); ++i)
        {
            std::cout << "Name: " << group_[i].get_name() << std::endl;
            std::cout << "Lastname: " << group_[i].get_lastname() << std::endl;
            std::cout << "Group: " << student.get_name_group() << std::endl;

            std::map<std::string, std::vector<uint16_t>> estimations_ = group_[i].get_estimations();

            for (auto it : estimations_)
            {
                std::cout << it.first << "\tEstimations: ";
                for (auto es : it.second)
                {
                    std::cout  << es << " ";
                }
            }
        }
	}
private:
    Studet student;
    std::string name_group_;
    std::vector<Studet> group_;
};
