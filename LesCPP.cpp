﻿#include <iostream>
#include <stack>

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
    std::auto_ptr <T*> current_array_;
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
        for (size_t i = 0; i < size; i++)
        {
            current_array_[i] = generate_random_number();
        }
    }
    Arr() = default;

    T find_min() const
    {
        if (size_ <= 0) return -1;
        T min = current_array_[0];
        for (size_t i = 1; i < size_; ++i) {
            if (current_array_[i] < min) {
                min = current_array_[i];
            }
        }
        return min;
    }

    T find_max() const
    {
        if (size_ <= 0) return -1;
        T max = current_array_[0];
        for (size_t i = 1; i < size_; ++i) {
            if (current_array_[i] > max) {
                max = current_array_[i];
            }
        }
        return max;
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

//Необходимо написать функцию, которая принимает входную строку, состоящую из скобок(, ), { , }, [, ], и проверяет, является ли данная строка сбалансированной.
//Строка считается сбалансированной, если каждая открывающая скобка имеет соответствующую закрывающую скобку в правильном порядке.

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

int main()
{
    std::string input_str = "((([{([)}])])(){({){}})}{{{(}}{[]})}";
    check_balance(input_str);
}
